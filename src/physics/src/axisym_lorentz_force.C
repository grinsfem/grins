//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "grins/axisym_lorentz_force.h"

// GRINS
#include "grins_config.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_context.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace GRINS
{

  AxisymmetricLorentzForce::AxisymmetricLorentzForce( const std::string& physics_name,
						      const GetPot& input )
    : Physics(physics_name,input),
      _factor( input("Physics/"+physics_name+"/factor", 1.0 ) )
  {
    this->read_input_options(input);
    return;
  }

  AxisymmetricLorentzForce::~AxisymmetricLorentzForce()
  {
    return;
  }

  void AxisymmetricLorentzForce::read_input_options( const GetPot& input )
  {
    this->_u_FE_family =
      libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+axisymmetric_incomp_navier_stokes+"/FE_family", "LAGRANGE") );

    this->_u_order =
      libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+axisymmetric_incomp_navier_stokes+"/u_order", "SECOND") );

    this->_V_FE_family =
      libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+axisymmetric_electrostatics+"/FE_family", "LAGRANGE") );

    this->_V_order =
      libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+axisymmetric_electrostatics+"/V_order", "FIRST") );

    this->_A_FE_family =
      libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+axisymmetric_magnetostatics+"/FE_family", "NEDELEC_ONE") );

    this->_A_order =
      libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+axisymmetric_magnetostatics+"/A_order", "FIRST") );

    // Read variable naming info
    this->_u_r_var_name = input("Physics/VariableNames/r_velocity", u_r_var_name_default );
    this->_u_z_var_name = input("Physics/VariableNames/z_velocity", u_z_var_name_default );
    this->_A_var_name = input("Physics/VariableNames/MagneticPotential", A_var_name_default );
    this->_V_var_name = input("Physics/VariableNames/ElectricPotential", V_var_name_default );

    _sigma = input("Physics/"+axisymmetric_magnetostatics+"/sigma", 1.0);

    return;
  }

  void AxisymmetricLorentzForce::init_variables( libMesh::FEMSystem* system )
  {
    _u_r_var = system->add_variable(_u_r_var_name, _u_order, _u_FE_family);
    _u_z_var = system->add_variable(_u_z_var_name, _u_order, _u_FE_family);
    _A_var   = system->add_variable(_A_var_name, _A_order, _A_FE_family);
    _V_var   = system->add_variable(_V_var_name, _V_order, _V_FE_family);
    return;
  }

  void AxisymmetricLorentzForce::element_time_derivative( bool request_jacobian,
							  libMesh::FEMContext& context )
  {
#ifdef USE_GRVY_TIMERS
    this->_timer->BeginTimer("AxisymmetricLorentzForce::element_time_derivative");
#endif

    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.dof_indices_var[_u_r_var].size();
    const unsigned int n_A_dofs = context.dof_indices_var[_A_var].size();
    const unsigned int n_V_dofs = context.dof_indices_var[_V_var].size();

    // Get finite element object
    FEGenericBase<RealGradient>* A_fe;
    FEGenericBase<Real>* V_fe;
    FEGenericBase<Real>* u_r_fe;
    context.get_element_fe<RealGradient>( _A_var, A_fe );
    context.get_element_fe<Real>( _V_var, V_fe );
    context.get_element_fe<Real>( _u_r_var, u_r_fe );

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      A_fe->get_JxW();

    // The velocity shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& A_curl_phi =
      A_fe->get_curl_phi();

    const std::vector<std::vector<libMesh::RealGradient> >& V_gradphi =
      V_fe->get_dphi();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& vel_phi =
      u_r_fe->get_phi();

    // Physical location of the quadrature points
    const std::vector<libMesh::Point>& qpoint =
      u_r_fe->get_xyz();

    // Get residuals
    libMesh::DenseSubVector<Number> &Fr = *context.elem_subresiduals[_u_r_var]; // R_{r}
    libMesh::DenseSubVector<Number> &Fz = *context.elem_subresiduals[_u_z_var]; // R_{z}

    // Get Jacobians
    libMesh::DenseSubMatrix<Number> &KrV = *context.elem_subjacobians[_u_r_var][_V_var]; // R_{r},{T}
    libMesh::DenseSubMatrix<Number> &KzV = *context.elem_subjacobians[_u_z_var][_V_var]; // R_{z},{T}
    libMesh::DenseSubMatrix<Number> &KrA = *context.elem_subjacobians[_u_r_var][_A_var]; // R_{r},{T}
    libMesh::DenseSubMatrix<Number> &KzA = *context.elem_subjacobians[_u_z_var][_A_var]; // R_{z},{T}

    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.element_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	const libMesh::Number r = qpoint[qp](0);

	// Compute the solution & its gradient at the old Newton iterate.
	libMesh::Gradient grad_V, curl_A;

	context.interior_gradient(_V_var, qp, grad_V);
	context.interior_curl(_A_var, qp, curl_A);

	// First, an i-loop over the velocity degrees of freedom.
	// We know that n_u_dofs == n_v_dofs so we can compute contributions
	// for both at the same time.
	for (unsigned int i=0; i != n_u_dofs; i++)
	  {
	    Fr(i) += _factor*_sigma*grad_V(1)*curl_A(2)*vel_phi[i][qp]*r*JxW[qp];
	    Fz(i) += -_factor*_sigma*grad_V(0)*curl_A(2)*vel_phi[i][qp]*r*JxW[qp];

	    if (request_jacobian)
	      {
		for (unsigned int j=0; j != n_V_dofs; j++)
		  {
		    KrV(i,j) += _factor*_sigma*V_gradphi[j][qp](1)*curl_A(2)*vel_phi[i][qp]*r*JxW[qp];
		    KzV(i,j) += -_factor*_sigma*V_gradphi[j][qp](0)*curl_A(2)*vel_phi[i][qp]*r*JxW[qp];
		  } // End j dof loop
	      
		for (unsigned int j=0; j != n_A_dofs; j++)
		  {
		    KrA(i,j) += _factor*_sigma*grad_V(1)*A_curl_phi[j][qp](2)*vel_phi[i][qp]*r*JxW[qp];
		    KzA(i,j) += -_factor*_sigma*grad_V(0)*A_curl_phi[j][qp](2)*vel_phi[i][qp]*r*JxW[qp];
		  } // End j dof loop

	      } // End request_jacobian check

	  } // End i dof loop
      } // End quadrature loop

#ifdef USE_GRVY_TIMERS
    this->_timer->EndTimer("AxisymmetricLorentzForce::element_time_derivative");
#endif

    return;
  }

} // namespace GRINS
