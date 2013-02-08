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
#include "grins_config.h"
#include "grins/boussinesq_buoyancy.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_context.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace GRINS
{

  BoussinesqBuoyancy::BoussinesqBuoyancy( const std::string& physics_name, const GetPot& input )
    : Physics(physics_name,input),
      _T_FE_family( libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+heat_transfer+"/FE_family", "LAGRANGE") ) ),
      _V_FE_family( libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+incompressible_navier_stokes+"/FE_family", "LAGRANGE") ) ),
      _T_order( libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+heat_transfer+"/T_order", "SECOND") ) ),
      _V_order( libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+incompressible_navier_stokes+"/V_order", "SECOND") ) ),
      _u_var_name( input("Physics/VariableNames/u_velocity", u_var_name_default ) ),
      _v_var_name( input("Physics/VariableNames/v_velocity", v_var_name_default ) ),
      _w_var_name( input("Physics/VariableNames/w_velocity", w_var_name_default ) ),
      _T_var_name( input("Physics/VariableNames/Temperature", T_var_name_default ) ),
      _rho_ref( input("Physics/"+boussinesq_buoyancy+"/rho_ref", 1.0) ),
      _T_ref( input("Physics/"+boussinesq_buoyancy+"/T_ref", 1.0) ),
      _beta_T( input("Physics/"+boussinesq_buoyancy+"/beta_T", 1.0) )
  {
    unsigned int g_dim = input.vector_variable_size("Physics/"+boussinesq_buoyancy+"/g");

    _g(0) = input("Physics/"+boussinesq_buoyancy+"/g", 0.0, 0 );
    _g(1) = input("Physics/"+boussinesq_buoyancy+"/g", 0.0, 1 );
  
    if( g_dim == 3)
      _g(2) = input("Physics/"+boussinesq_buoyancy+"/g", 0.0, 2 );

    return;
  }

  BoussinesqBuoyancy::~BoussinesqBuoyancy()
  {
    return;
  }

  void BoussinesqBuoyancy::init_variables( libMesh::FEMSystem* system )
  {
    // Get libMesh to assign an index for each variable
    this->_dim = system->get_mesh().mesh_dimension();

    _T_var = system->add_variable( _T_var_name, this->_T_order, _T_FE_family);
 
    // If these are already added, then we just get the index. 
    _u_var = system->add_variable(_u_var_name, _V_order, _V_FE_family );
    _v_var = system->add_variable(_v_var_name, _V_order, _V_FE_family );
    if (_dim == 3)
      _w_var = system->add_variable(_w_var_name, _V_order, _V_FE_family );

    return;
  }

  void BoussinesqBuoyancy::element_time_derivative( bool compute_jacobian,
						    libMesh::FEMContext& context,
						    CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("BoussinesqBuoyancy::element_time_derivative");
#endif
  
    if (_dim != 3)
      _w_var = _u_var; // for convenience

    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.dof_indices_var[_u_var].size();
    const unsigned int n_T_dofs = context.dof_indices_var[_T_var].size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.element_fe_var[_u_var]->get_JxW();

    // The velocity shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& vel_phi =
      context.element_fe_var[_u_var]->get_phi();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.element_fe_var[_T_var]->get_phi();

    // Get residuals
    libMesh::DenseSubVector<libMesh::Number> &Fu = *context.elem_subresiduals[_u_var]; // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = *context.elem_subresiduals[_v_var]; // R_{v}
    libMesh::DenseSubVector<libMesh::Number> &Fw = *context.elem_subresiduals[_w_var]; // R_{w}

    // Get Jacobians
    libMesh::DenseSubMatrix<libMesh::Number> &KuT = *context.elem_subjacobians[_u_var][_T_var]; // R_{u},{T}
    libMesh::DenseSubMatrix<libMesh::Number> &KvT = *context.elem_subjacobians[_v_var][_T_var]; // R_{v},{T}
    libMesh::DenseSubMatrix<libMesh::Number> &KwT = *context.elem_subjacobians[_w_var][_T_var]; // R_{w},{T}

    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.element_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	// Compute the solution & its gradient at the old Newton iterate.
	libMesh::Number T;
	T = context.interior_value(_T_var, qp);

	// First, an i-loop over the velocity degrees of freedom.
	// We know that n_u_dofs == n_v_dofs so we can compute contributions
	// for both at the same time.
	for (unsigned int i=0; i != n_u_dofs; i++)
	  {
	    Fu(i) += -_rho_ref*_beta_T*(T - _T_ref)*_g(0)*vel_phi[i][qp]*JxW[qp];
	    Fv(i) += -_rho_ref*_beta_T*(T - _T_ref)*_g(1)*vel_phi[i][qp]*JxW[qp];

	    if (_dim == 3)
	      Fw(i) += -_rho_ref*_beta_T*(T - _T_ref)*_g(2)*vel_phi[i][qp]*JxW[qp];

	    if (compute_jacobian)
	      {
		for (unsigned int j=0; j != n_T_dofs; j++)
		  {
		    KuT(i,j) += -_rho_ref*_beta_T*_g(0)*vel_phi[i][qp]*T_phi[j][qp]*JxW[qp];
		    KvT(i,j) += -_rho_ref*_beta_T*_g(1)*vel_phi[i][qp]*T_phi[j][qp]*JxW[qp];

		    if (_dim == 3)
		      KwT(i,j) += -_rho_ref*_beta_T*_g(2)*vel_phi[i][qp]*T_phi[j][qp]*JxW[qp];

		  } // End j dof loop
	      } // End compute_jacobian check

	  } // End i dof loop
      } // End quadrature loop

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("BoussinesqBuoyancy::element_time_derivative");
#endif

    return;
  }

} // namespace GRINS
