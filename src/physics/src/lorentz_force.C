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

// This class
#include "grins/lorentz_force.h"
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace GRINS
{
  LorentzForce::LorentzForce( const std::string& physics_name, const GetPot& input )
    : Physics(physics_name,input),
      _u_FE_family( libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+incompressible_navier_stokes+"/FE_family", "LAGRANGE") ) ),
      _A_FE_family( libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+magnetostatics+"/FE_family", "NEDELEC_ONE") ) ),
      _V_FE_family( libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+electrostatics+"/FE_family", "LAGRANGE") ) ),
      _u_order( libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+incompressible_navier_stokes+"/u_order", "SECOND") ) ),
      _A_order( libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+magnetostatics+"/A_order", "FIRST") ) ),
      _V_order( libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+electrostatics+"/V_order", "FIRST") ) ),
      _u_var_name( input("Physics/VariableNames/u_velocity", u_var_name_default ) ),
      _v_var_name( input("Physics/VariableNames/v_velocity", v_var_name_default ) ),
      _w_var_name( input("Physics/VariableNames/w_velocity", w_var_name_default ) ),
      _A_var_name( input("Physics/VariableNames/MagneticPotential", A_var_name_default ) ),
      _V_var_name( input("Physics/VariableNames/ElectricPotential", V_var_name_default ) ),
      _sigma( input("Physics/"+magnetostatics+"/sigma", 1.0) )
  {
    return;
  }

  LorentzForce::~LorentzForce()
  {
    return;
  }

  void LorentzForce::init_variables( libMesh::FEMSystem* system )
  {
    // Get libMesh to assign an index for each variable
    this->_dim = system->get_mesh().mesh_dimension();

    _u_var = system->add_variable( _u_var_name, this->_u_order, _u_FE_family);
    _v_var = system->add_variable( _v_var_name, this->_u_order, _u_FE_family);

    if (_dim == 3)
      _w_var = system->add_variable( _w_var_name, this->_u_order, _u_FE_family);

    _A_var   = system->add_variable(_A_var_name, _A_order, _A_FE_family);
    _V_var   = system->add_variable(_V_var_name, _V_order, _V_FE_family);

    return;
  }

  void LorentzForce::element_time_derivative( bool compute_jacobian,
					      AssemblyContext& context,
					      CachedValues& cache )
  {
#ifdef USE_GRVY_TIMERS
    this->_timer->BeginTimer("LorentzForce::element_time_derivative");
#endif

    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.get_dof_indices(_u_var).size();
    const unsigned int n_A_dofs = context.get_dof_indices(_A_var).size();
    const unsigned int n_V_dofs = context.get_dof_indices(_V_var).size();

    // Get finite element object
    libMesh::FEGenericBase<libMesh::RealGradient>* A_fe;
    libMesh::FEGenericBase<libMesh::Real>* V_fe;
    libMesh::FEGenericBase<libMesh::Real>* u_fe;
    context.get_element_fe<libMesh::RealGradient>( _A_var, A_fe );
    context.get_element_fe<libMesh::Real>( _V_var, V_fe );
    context.get_element_fe<libMesh::Real>( _u_var, u_fe );

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real>& JxW = A_fe->get_JxW();

    // The velocity shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& A_curl_phi = A_fe->get_curl_phi();

    const std::vector<std::vector<libMesh::RealGradient> >& V_gradphi = V_fe->get_dphi();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& u_phi = u_fe->get_phi();

    // Get residuals
    libMesh::DenseSubVector<libMesh::Number>& Fu = context.get_elem_residual(_u_var); // R_{u}
    libMesh::DenseSubVector<libMesh::Number>& Fv = context.get_elem_residual(_v_var); // R_{v}
    libMesh::DenseSubVector<libMesh::Number>* Fw = NULL;

    // Get Jacobians
    libMesh::DenseSubMatrix<libMesh::Number>& KuV = context.get_elem_jacobian(_u_var,_V_var); // R_{u},{V}
    libMesh::DenseSubMatrix<libMesh::Number>& KvV = context.get_elem_jacobian(_v_var,_V_var); // R_{v},{V}
    libMesh::DenseSubMatrix<libMesh::Number>* KwV = NULL;

    libMesh::DenseSubMatrix<libMesh::Number>& KuA = context.get_elem_jacobian(_u_var,_A_var); // R_{u},{A}
    libMesh::DenseSubMatrix<libMesh::Number>& KvA = context.get_elem_jacobian(_v_var,_A_var); // R_{v},{A}
    libMesh::DenseSubMatrix<libMesh::Number>* KwA = NULL;

    if( this->_dim == 3)
      {
        Fw  = &context.get_elem_residual(_w_var); // R_{w}
        KwV = &context.get_elem_jacobian(_w_var,_V_var); // R_{w},{V}
        KwA = &context.get_elem_jacobian(_w_var,_A_var); // R_{w},{A}
      }

    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	// Compute the solution & its gradient at the old Newton iterate.
	libMesh::Gradient grad_V, curl_A;

	context.interior_gradient(_V_var, qp, grad_V);
	context.interior_curl(_A_var, qp, curl_A);

	libMesh::Gradient lorentz_force = -_sigma*grad_V.cross(curl_A);

	// First, an i-loop over the velocity degrees of freedom.
	// We know that n_u_dofs == n_v_dofs so we can compute contributions
	// for both at the same time.
	for (unsigned int i=0; i != n_u_dofs; i++)
	  {
	    Fu(i) += lorentz_force(0)*u_phi[i][qp]*JxW[qp];
	    Fv(i) += lorentz_force(1)*u_phi[i][qp]*JxW[qp];

	    if(_dim == 3 )
	      {
		(*Fw)(i) += lorentz_force(2)*u_phi[i][qp]*JxW[qp];
	      }

	    if (compute_jacobian)
	      {
		for (unsigned int j=0; j != n_V_dofs; j++)
		  {
		    KuV(i,j) += u_phi[i][qp]*JxW[qp]*(-_sigma)*( V_gradphi[j][qp](1)*curl_A(2) 
								 - V_gradphi[j][qp](2)*curl_A(1) );
		    
		    KvV(i,j) += u_phi[i][qp]*JxW[qp]*(_sigma)*( V_gradphi[j][qp](0)*curl_A(2) 
								- V_gradphi[j][qp](2)*curl_A(0) );
		    if(_dim == 3 )
		      {
			(*KwV)(i,j) += u_phi[i][qp]*JxW[qp]*(-_sigma)*( V_gradphi[j][qp](0)*curl_A(1) 
								     - V_gradphi[j][qp](1)*curl_A(0) );
		      }
		  } // End j dof loop
	      
		for (unsigned int j=0; j != n_A_dofs; j++)
		  {
		    KuA(i,j) += u_phi[i][qp]*JxW[qp]*(-_sigma)*( grad_V(1)*A_curl_phi[j][qp](2) 
								 - grad_V(2)*A_curl_phi[j][qp](1) );

		    KvA(i,j) += u_phi[i][qp]*JxW[qp]*(_sigma)*( grad_V(0)*A_curl_phi[j][qp](2) 
								- grad_V(2)*A_curl_phi[j][qp](0) );

		    if(_dim == 3 )
		      {
			(*KwA)(i,j) += u_phi[i][qp]*JxW[qp]*(-_sigma)*( grad_V(0)*A_curl_phi[j][qp](1) 
								     - grad_V(1)*A_curl_phi[j][qp](0) );
		      }
		  } // End j dof loop

	      } // End request_jacobian check

	  } // End i dof loop
      } // End quadrature loop

#ifdef USE_GRVY_TIMERS
    this->_timer->EndTimer("LorentzForce::element_time_derivative");
#endif

    return;
  }

  void LorentzForce::compute_element_cache( const AssemblyContext& context, 
					    const std::vector<libMesh::Point>& points,
					    CachedValues& cache )
  {
    // Lorentz Force
    if( cache.is_active(Cache::LORENTZ_FORCE_X) )
      {
	std::vector<libMesh::Real> Lx, Ly, Lz;
	Lx.reserve( points.size() );
	Ly.reserve( points.size() );
	if( _dim > 2 )
	  Lz.reserve( points.size() );

	for( std::vector<libMesh::Point>::const_iterator point = points.begin();
	     point != points.end(); point++ )
	  {
	    libMesh::Gradient J, B;
	    context.point_curl(_A_var,* point, B);
	    context.point_gradient(_V_var, *point, J);
	    // J = \sigma*E = -\sigma \nabla V
	    J *= -_sigma;
	    libMesh::Gradient L = J.cross(B);

	    Lx.push_back(L(0));
	    Ly.push_back(L(1));
	    if( _dim > 2 )
	      Lz.push_back(L(2));
	  }

	cache.set_values( Cache::LORENTZ_FORCE_X, Lx );
	cache.set_values( Cache::LORENTZ_FORCE_Y, Ly );
	if( _dim > 2 )
	  cache.set_values( Cache::LORENTZ_FORCE_Z, Lz );
      }

    return;
  }

} // namespace GRINS
