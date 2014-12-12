//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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
#include "grins/spalart_allmaras.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/generic_ic_handler.h"
#include "grins/inc_navier_stokes_bc_handling.h"

#include "grins/constant_viscosity.h"
#include "grins/parsed_viscosity.h"
#include "grins/turbulence_viscosity.h" 

#include "grins/turbulence_models_macro.h"

// libMesh
#include "libmesh/quadrature.h"
#include "libmesh/elem.h"

namespace GRINS
{

  template<class Mu>
  SpalartAllmaras<Mu>::SpalartAllmaras(const std::string& physics_name, const GetPot& input )
    : Physics(physics_name, input), // Define class variables
      _turbulence_vars(input, spalartallmaras),
      _cb1(0.1355),
      _sigma(2./3.),
      _cb2(0.622),
      _kappa(0.41)
  {    
    // This is deleted in the base class
    this->_bc_handler = new IncompressibleNavierStokesBCHandling( physics_name, input );

    if( this->_bc_handler->is_axisymmetric() )
      {
        this->_is_axisymmetric = true;
      }

    this->_ic_handler = new GenericICHandler( physics_name, input );

    return;
  }

  template<class Mu>
  SpalartAllmaras<Mu>::~SpalartAllmaras()
  {
    return;
  }

  template<class Mu>
  void SpalartAllmaras<Mu>::init_variables( libMesh::FEMSystem* system )
  {
    this->_dim = system->get_mesh().mesh_dimension();
    
    this->_turbulence_vars.init(system); // Should replace this turbulence_vars

    return;
  }

  template<class Mu>
  void SpalartAllmaras<Mu>::init_context( AssemblyContext& context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    context.get_element_fe(_turbulence_vars._nu_var())->get_JxW();
    context.get_element_fe(_turbulence_vars._nu_var())->get_phi();
    context.get_element_fe(_turbulence_vars._nu_var())->get_dphi();
    context.get_element_fe(_turbulence_vars._nu_var())->get_xyz();

    context.get_element_fe(_turbulence_vars._nu_var())->get_phi();
    context.get_element_fe(_turbulence_vars._nu_var())->get_xyz();

    context.get_side_fe(_turbulence_vars._nu_var())->get_JxW();
    context.get_side_fe(_turbulence_vars._nu_var())->get_phi();
    context.get_side_fe(_turbulence_vars._nu_var())->get_dphi();
    context.get_side_fe(_turbulence_vars._nu_var())->get_xyz();

    return;
  }

  template<class Mu>
  void SpalartAllmaras<Mu>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    const unsigned int dim = system->get_mesh().mesh_dimension();

    // Tell the system to march velocity forward in time, but
    // leave p as a constraint only
    system->time_evolving(_turbulene_vars._nu_var());
    
    return;
  }
 
  template<class Mu>
  void SpalartAllmaras<Mu>::element_time_derivative( bool compute_jacobian,
                                                            AssemblyContext& context,
                                                            CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("SpalartAllmaras::element_time_derivative");
#endif

    // Get a pointer to the current element, we need this for computing the distance to wall for the
    // quadrature points
    Elem &elem_pointer = context.get_elem();

    // The number of local degrees of freedom in each variable.
    const unsigned int n_nu_dofs = context.get_dof_indices(this->_turbulence_vars.nu_var()).size();
        
    // We get some references to cell-specific data that
    // will be used to assemble the linear system.

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_turbulence_vars.nu_var())->get_JxW();

    // The viscosity shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& nu_phi =
      context.get_element_fe(this->_turbulence_vars.nu_var())->get_phi();

    // The viscosity shape function gradients (in global coords.)
    // at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& nu_gradphi =
      context.get_element_fe(this->_turbulence_vars.nu_var())->get_dphi();
    
    // Quadrature point locations
    const std::vector<libMesh::Point>& nu_qpoint = 
      context.get_element_fe(this->_turbulence_vars.nu_var())->get_xyz();

    // The subvectors and submatrices we need to fill:
    //
    // K_{\alpha \beta} = R_{\alpha},{\beta} = \partial{ R_{\alpha} } / \partial{ {\beta} } (where R denotes residual)
    // e.g., for \alpha = v and \beta = u we get: K{vu} = R_{v},{u}
    // Note that Kpu, Kpv, Kpw and Fp comes as constraint.

    libMesh::DenseSubMatrix<libMesh::Number> &Knunu = context.get_elem_jacobian(this->_turbulence_vars.nu_var(), this->_turbulence_vars.nu_var()); // R_{nu},{nu}
    
    libMesh::DenseSubVector<libMesh::Number> &Fnu = context.get_elem_residual(this->_turbulence_vars.nu_var()); // R_{nu}
    
    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.get_element_qrule().n_points();

    // Auto pointer to distance fcn evaluated at quad points
    AutoPtr< DenseVector<Real> > distance_qp;

    // Fill the vector of distances to quadrature points
    distance_qp = this->distance_function->interpolate(elem_pointer, context.get_element_qrule().get_points());

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // Compute the solution & its gradient at the old Newton iterate.
        libMesh::Number nu;
        nu = context.interior_value(this->_turbulence_vars.nu_var(), qp);        

        libMesh::Gradient grad_nu;
        grad_nu = context.interior_gradient(this->_turbulence_vars.nu_var(), qp);
        
        const libMesh::Number  grad_nu_x = grad_nu(0);
        const libMesh::Number  grad_nu_y = grad_nu(1);
        const libMesh::Number  grad_nu_z = (this->_dim == 3)?grad_nu(2):0;
        
        const libMesh::Number x = u_qpoint[qp](0);
	const libMesh::Number y = u_qpoint[qp](1);
	const libMesh::Number z = (this->_dim==3)?u_qpoint[qp](2):0;

        libMesh::Real jac = JxW[qp];
	
	// The physical viscosity
	libMesh::Real _mu_qp = this->_mu(context, qp);
   
        // First, an i-loop over the viscosity degrees of freedom.        
        for (unsigned int i=0; i != n_nu_dofs; i++)
          {
	    // TODO: intialize constants cb1, cb2, cw1, sigma, and functions source_fn(nu), destruction_fn(nu), and resolve issue of grad(nu + nu_tilde)
            Fnu(i) += jac *
              (-this->_cb1*this->_source_fn(nu, _mu_qp, distance_qp[qp])*nu*nu_phi[i][qp]  // source term
	       + (1./_sigma)*(-(_mu_qp+nu)*grad_nu*nu_gradphi[i][qp] - grad_nu*grad_nu*phi[i][qp] + this->_cb2*grad_nu*grad_nu*phi[i][qp])  // diffusion term 
               - this->_cw1*this->_destruction_fn(nu)*pow(nu/x, 2.)*phi[i][qp]) // destruction term                          
	      // Compute the jacobian if not using numerical jacobians  
	      if (compute_jacobian)
		{
		  for (unsigned int j=0; j != n_u_dofs; j++)
		    {
		      // TODO: precompute some terms like:
		      //   (Uvec*u_gradphi[j][qp]),
		      //   u_phi[i][qp]*u_phi[j][qp],
		      //   (u_gradphi[i][qp]*u_gradphi[j][qp])
		      
		      Knunu(i,j) += jac *
			(      // source term
			 // diffusion term
			       ); // destruction term
                                                              
		    } // end of the inner dof (j) loop                                

		} // end - if (compute_jacobian)

          } // end of the outer dof (i) loop
      } // end of the quadrature point (qp) loop

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("SpalartAllmaras::element_time_derivative");
#endif

    return;
  }
    
  template<class Mu>
  Real SpalartAllmaras<Mu>::_source_fn(libMesh::Number nu, libMesh::Real mu, Real wall_distance)
  {
    // Step 1
    Real _kai = nu/mu;
    
    // Step 2
    Real _fv1 = pow(_kai, 3.0)/(pow(_kai, 3.0) + pow(_cv1, 3.0));

    // Step 3
    Real _fv2 = 1 - (_kai/(1 + _kai*_fv1));

    // Step 4
    Real _S_bar = nu/(pow(_kappa, 2.0) * pow(wall_distance, 2.0)) ;

    // Step 5, the absolute value of the vorticity
    Real _S = fabs(_vorticity(FlowVars, qp));

    // Step 6
    Real _S_tilde = 0.0;
    if(_S_bar >= -_cv2*_S)
      {
	_S_tilde = _S + _S_bar;
      }
    else
      {
	_S_tilde = _S + (_S*(pow(_cv2,2.0)*_S + _cv3*_S_bar))/((_cv3 - (2*_cv2))*_S - _S_bar);
      }
    
    return _S_tilde;
  }

  template<class Mu>
  Real SpalartAllmaras<Mu>::_destruction_fn(libMesh::Number nu, Real wall_distance, Real _S_tilde)
  {
    // Step 1
    Real _r = 0.0;
    
    if(nu/(_S_tilde*pow(_kappa,2.0)*pow(wall_distance,2.0)) < _r_lin)
      {
	_r = nu/(_S_tilde*pow(_kappa,2.0)*pow(wall_distance,2.0));
      }
    else
      {
	_r = _r_lin;
      }

    // Step 2
    Real _g = _r + _c_w2*(pow(_r,6.0) - _r);

    // Step 3
    Real _fw = _g*pow((1 + pow(_c_w3,6.0))/(pow(_g,6.0) + pow(_c_w3,6.0)), 1.0/6.0);
    
    return _fw;
  }
  

  
} // namespace GRINS

// Instantiate
INSTANTIATE_TURBULENCE_MODELS_SUBCLASS(SpalartAllmaras);
