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
#include "grins/magnetostatics.h"

// GRINS
#include "grins/magnetostatics_bc_handling.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_context.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace GRINS
{
  Magnetostatics::Magnetostatics( const std::string& physics_name, const GetPot& input )
    : Physics(physics_name,input),
      _A_var_name( input("Physics/VariableNames/MagneticPotential", A_var_name_default ) ),
      _V_var_name( input("Physics/VariableNames/ElectricPotential", V_var_name_default ) ),
      _A_FE_family( libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+magnetostatics+"/FE_family", "NEDELEC_ONE") ) ),
      _V_FE_family( libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+electrostatics+"/FE_family", "LAGRANGE") ) ),
      _A_order( libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+magnetostatics+"/A_order", "FIRST") ) ),
      _V_order( libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+electrostatics+"/V_order", "SECOND") ) ),
      _sigma( input("Physics/"+magnetostatics+"/sigma", 1.0 ) ),
      _mu( input("Physics/"+magnetostatics+"/mu", 1.0 ) )
  {
    // This is deleted in the base class
    _bc_handler = new GRINS::MagnetostaticsBCHandling( physics_name, input );

    return;
  }

  Magnetostatics::~Magnetostatics()
  {
    return;
  }

  void Magnetostatics::init_variables( libMesh::FEMSystem* system )
  {
    // Get libMesh to assign an index for each variable
    this->_dim = system->get_mesh().mesh_dimension();
    
    _A_var = system->add_variable( _A_var_name, _A_order, _A_FE_family);
    
    _V_var = system->add_variable( _V_var_name, _V_order, _V_FE_family);

    return;
  }

  void Magnetostatics::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    system->time_evolving(_A_var);

    return;
  }

  void Magnetostatics::init_context( libMesh::FEMContext& context )
  {   
    libMesh::FEGenericBase<libMesh::RealGradient>* A_fe;
    libMesh::FEGenericBase<libMesh::Real>* V_fe;
    context.get_element_fe<libMesh::RealGradient>( _A_var, A_fe );
    context.get_element_fe<libMesh::Real>( _V_var, V_fe );

    libMesh::FEGenericBase<libMesh::RealGradient>* A_side_fe;
    libMesh::FEGenericBase<libMesh::Real>* V_side_fe;
    context.get_side_fe<libMesh::RealGradient>( _A_var, A_side_fe );
    context.get_side_fe<libMesh::Real>( _V_var, V_side_fe );

    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    A_fe->get_JxW();
    A_fe->get_phi();
    A_fe->get_curl_phi();
    A_fe->get_xyz();
    
    V_fe->get_JxW();
    V_fe->get_phi();
    V_fe->get_dphi();
    V_fe->get_xyz();
    
    A_side_fe->get_JxW();
    A_side_fe->get_phi();
    A_side_fe->get_curl_phi();
    A_side_fe->get_xyz();
    
    V_side_fe->get_JxW();
    V_side_fe->get_phi();
    V_side_fe->get_dphi();
    V_side_fe->get_xyz();
    
    return;
  }

  void Magnetostatics::element_time_derivative( bool compute_jacobian,
						libMesh::FEMContext& context,
						CachedValues& cache )
  {
#ifdef USE_GRVY_TIMERS
    this->_timer->BeginTimer("Magnetostatics::element_time_derivative");
#endif
    
    // The number of local degrees of freedom in each variable.
    const unsigned int n_A_dofs = context.dof_indices_var[_A_var].size();
    const unsigned int n_V_dofs = context.dof_indices_var[_V_var].size();
  
    // Get finite element object
    libMesh::FEGenericBase<libMesh::RealGradient>* A_fe;
    libMesh::FEGenericBase<libMesh::Real>* V_fe;
    context.get_element_fe<libMesh::RealGradient>( _A_var, A_fe );
    context.get_element_fe<libMesh::Real>( _V_var, V_fe );

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real>& JxW = A_fe->get_JxW();

    // The magnetic potential shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& A_phi = A_fe->get_phi();
    
    // The electric potential shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& A_curl_phi = A_fe->get_curl_phi();
    
    // The electric potential shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& V_gradphi = V_fe->get_dphi();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<libMesh::Number> &F_A = *context.elem_subresiduals[_A_var]; // R_{A}
    
    libMesh::DenseSubMatrix<libMesh::Number> &K_AA = *context.elem_subjacobians[_A_var][_A_var]; // R_{A},{A}

    libMesh::DenseSubMatrix<libMesh::Number> &K_AV = *context.elem_subjacobians[_A_var][_V_var]; // R_{A},{V}

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
	libMesh::Gradient curl_A, grad_V;

	context.interior_gradient(_V_var, qp, grad_V);
	context.interior_curl(_A_var, qp, curl_A);

	// Loop over magnetic potential dofs
	for (unsigned int i=0; i != n_A_dofs; i++)
	  {
	    F_A(i) += (curl_A*A_curl_phi[i][qp]/_mu + _sigma*grad_V*A_phi[i][qp])*JxW[qp];

	    if( compute_jacobian )
	      {
		for (unsigned int j=0; j != n_V_dofs; j++)
		  {
		    K_AV(i,j) += _sigma*V_gradphi[j][qp]*A_phi[i][qp]*JxW[qp];
		  }

		for (unsigned int j=0; j != n_A_dofs; j++)
		  {
		    K_AA(i,j) += A_curl_phi[j][qp]*A_curl_phi[i][qp]/_mu*JxW[qp];
		  }
	      }
	  } // magnetic potential dofs

      } // end quadrature loop

    return;
  }

  void Magnetostatics::side_time_derivative( bool compute_jacobian,
					     libMesh::FEMContext& context,
					     CachedValues& cache )
  {
#ifdef USE_GRVY_TIMERS
    this->_timer->BeginTimer("Magnetostatics::side_time_derivative");
#endif

    std::vector<BoundaryID> ids = context.side_boundary_ids();

    for( std::vector<BoundaryID>::const_iterator it = ids.begin();
         it != ids.end(); it++ )
      {
        libmesh_assert (*it != libMesh::BoundaryInfo::invalid_id);
    
	_bc_handler->apply_neumann_bcs( context, cache, compute_jacobian, *it );
      }

#ifdef USE_GRVY_TIMERS
    this->_timer->EndTimer("Magnetostatics::side_time_derivative");
#endif
    
    return;
  }
  
  void Magnetostatics::side_constraint( bool compute_jacobian,
					libMesh::FEMContext& context,
					CachedValues& cache )
  {
#ifdef USE_GRVY_TIMERS
    //this->_timer->BeginTimer("Magnetostatics::side_constraint");
#endif
    
    // Get finite element object
    libMesh::FEGenericBase<libMesh::RealGradient>* side_fe = NULL;
    context.get_side_fe<libMesh::RealGradient>( _A_var, side_fe );

    // First we get some references to cell-specific data that
    // will be used to assemble the linear system.
    
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<libMesh::Real> &JxW = side_fe->get_JxW();
    
    // The velocity shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& phi = side_fe->get_phi();
    
    // The number of local degrees of freedom in each variable
    const unsigned int n_A_dofs = context.dof_indices_var[_A_var].size();
    
    const std::vector<libMesh::Point>& normals = side_fe->get_normals();
    
    // The penalty value.  \frac{1}{\epsilon}
    // in the discussion above.
    const Real penalty = 1.e10;
    
    libMesh::DenseSubMatrix<libMesh::Number> &K = *context.elem_subjacobians[_A_var][_A_var];
    libMesh::DenseSubVector<libMesh::Number> &F = *context.elem_subresiduals[_A_var];
    
    const unsigned int n_qpoints = context.get_side_qrule().n_points();
    
    std::vector<BoundaryID> ids = context.side_boundary_ids();
    for( std::vector<BoundaryID>::const_iterator it = ids.begin();
         it != ids.end(); it++ )
      {
	for (unsigned int qp=0; qp != n_qpoints; qp++)
	  {
	    libMesh::Gradient A;
	    context.side_value<libMesh::RealGradient>( _A_var, qp, A );
	
	    /*
	    if( _bc_handler->get_dirichlet_bc_type(*it) == MagnetostaticsBCHandling::AXISYMMETRIC )
	      {
		for (unsigned int i=0; i != n_A_dofs; i++)
		  {
		    F(i) += penalty*(A*phi[i][qp])*JxW[qp];
		
		    if (request_jacobian)
		      {
			for (unsigned int j=0; j != n_A_dofs; j++)
			  K(i,j) += penalty*(phi[j][qp]*phi[i][qp])*JxW[qp];
		      }
		  }
	      }
	    else
	      {
		libMesh::RealGradient N( normals[qp](0), normals[qp](1) );
	    
		libMesh::Gradient NcA = A.cross(N);
	    
		for (unsigned int i=0; i != n_A_dofs; i++)
		  {
		    F(i) += penalty*NcA*(phi[i][qp].cross(N))*JxW[qp];
		
		    if (request_jacobian)
		      {
			for (unsigned int j=0; j != n_A_dofs; j++)
			  K(i,j) += penalty*(phi[j][qp].cross(N))*(phi[i][qp].cross(N))*JxW[qp];
		      }
		  }
	      }
	    */
	  } // End quadrature loop
      }

#ifdef USE_GRVY_TIMERS
    this->_timer->EndTimer("Magnetostatics::side_constraint");
#endif

    return;
  }
  

  void Magnetostatics::mass_residual( bool compute_jacobian,
				      libMesh::FEMContext& context,
				      CachedValues& cache )
  {
#ifdef USE_GRVY_TIMERS
    this->_timer->BeginTimer("Magnetostatics::mass_residual");
#endif
    libmesh_not_implemented();

    // Get finite element object
    libMesh::FEGenericBase<libMesh::RealGradient>* A_fe;
    libMesh::FEGenericBase<libMesh::Real>* V_fe;
    context.get_element_fe<libMesh::RealGradient>( _A_var, A_fe );
    context.get_element_fe<libMesh::Real>( _V_var, V_fe );

    // First we get some references to cell-specific data that
    // will be used to assemble the linear system.

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW = A_fe->get_JxW();
  
    // The magnetic potential shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& A_phi = A_fe->get_phi();
  
    // The electric potential shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& V_phi = V_fe->get_phi();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<Real> &F_A = *context.elem_subresiduals[_A_var];

    libMesh::DenseSubVector<Real> &F_V = *context.elem_subresiduals[_V_var];

    libMesh::DenseSubMatrix<Real> &M_AA = *context.elem_subjacobians[_A_var][_A_var];

    libMesh::DenseSubMatrix<Real> &M_VA = *context.elem_subjacobians[_V_var][_A_var];

    unsigned int n_qpoints = context.element_qrule->n_points();

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
	// For the mass residual, we need to be a little careful.
	// The time integrator is handling the time-discretization
	// for us so we need to supply M(u_fixed)*u for the residual.
	// u_fixed will be given by the fixed_interior_* functions
	// while u will be given by the interior_* functions.
	libMesh::Gradient A_dot;
	context.fixed_interior_value<libMesh::RealGradient>(_A_var, qp, A_dot);

	/*
	  for (unsigned int i = 0; i != n_A_dofs; ++i)
	  {

          if( request_jacobian )
	  {
	  for (unsigned int j=0; j != n_T_dofs; j++)
	  {
                    
	  }
	  }// End of check on Jacobian
          
	  } // End of element dof loop
	*/      

      } // End of the quadrature point loop

#ifdef USE_GRVY_TIMERS
    this->_timer->EndTimer("Magnetostatics::mass_residual");
#endif

    return;
  }

  void Magnetostatics::compute_element_cache( const libMesh::FEMContext& context, 
					      const std::vector<libMesh::Point>& points,
					      CachedValues& cache )
  {
    // Magnetic Field
    if( cache.is_active(Cache::MAGNETIC_FIELD_X) )
      {
	std::vector<libMesh::Real> Hx, Hy, Hz;
	Hx.reserve( points.size() );
	
	if( _dim > 2 )
	  {
	    Hy.reserve( points.size() );
	    Hz.reserve( points.size() );
	  }

	for( std::vector<libMesh::Point>::const_iterator point = points.begin();
	     point != points.end(); point++ )
	  {
	    libMesh::Gradient H;
	    context.point_curl(_A_var,*point, H);

	    // In 2-D, the magnetic field is only the out-of-plane component
	    if( _dim == 2 )
	      Hx.push_back(H(2)/_mu);

	    if( _dim > 2 )
	      {
		Hx.push_back(H(0)/_mu);
		Hy.push_back(H(1)/_mu);
		Hz.push_back(H(2)/_mu);
	      }
	  }

	cache.set_values( Cache::MAGNETIC_FIELD_X, Hx );
	
	if( _dim > 2 )
	  {
	    cache.set_values( Cache::MAGNETIC_FIELD_Y, Hy );
	    cache.set_values( Cache::MAGNETIC_FIELD_Z, Hz );
	  }
	  
      }

    // Magnetic Flux
    if( cache.is_active(Cache::MAGNETIC_FLUX_X) )
      {
	std::vector<libMesh::Real> Bx, By, Bz;
	Bx.reserve( points.size() );
	
	if( _dim > 2 )
	  {
	    By.reserve( points.size() );
	    Bz.reserve( points.size() );
	  }

	if( cache.is_active(Cache::MAGNETIC_FIELD_X) )
	  {
	    const std::vector<libMesh::Number>& Hx = 
	      cache.get_cached_values( Cache::MAGNETIC_FIELD_X );

	    const std::vector<libMesh::Number>* Hy;
	    const std::vector<libMesh::Number>* Hz;
	    if( _dim > 2 )
	      {
		Hy = &cache.get_cached_values( Cache::MAGNETIC_FIELD_Y );
		Hz = &cache.get_cached_values( Cache::MAGNETIC_FIELD_Z );
	      }

	    for( unsigned int p = 0; p < points.size(); p++ )
	      {
		// In 2-D, the magnetic field is only the out-of-plane component
		Bx.push_back(_mu*Hx[p]);
		
		if( _dim > 2 )
		  {
		    By.push_back(_mu*(*Hy)[p]);
		    Bz.push_back(_mu*(*Hz)[p]);
		  }
	      }
	  }
	else
	  {
	    for( std::vector<libMesh::Point>::const_iterator point = points.begin();
		 point != points.end(); point++ )
	      {
		libMesh::Gradient B;
		context.point_curl(_A_var, *point, B);
		
		// In 2-D, the magnetic field is only the out-of-plane component
		if( _dim == 2 )
		  Bx.push_back(B(2));
		
		if( _dim > 2 )
		  {
		    Bx.push_back(B(0));
		    By.push_back(B(1));
		    Bz.push_back(B(2));
		  }
	      }
	  }

	cache.set_values( Cache::MAGNETIC_FLUX_X, Bx );
	
	if( _dim > 2 )
	  {
	    cache.set_values( Cache::MAGNETIC_FLUX_Y, By );
	    cache.set_values( Cache::MAGNETIC_FLUX_Z, Bz );
	  }
      }

    return;
  }

} // namespace GRINS
