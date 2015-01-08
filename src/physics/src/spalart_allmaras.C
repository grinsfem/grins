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
#include "grins/turbulence_models_macro.h"

#include "grins/constant_viscosity.h"
#include "grins/parsed_viscosity.h"
#include "grins/spalart_allmaras_viscosity.h"

// libMesh
#include "libmesh/quadrature.h"
#include "libmesh/elem.h"
#include "libmesh/unstructured_mesh.h"

namespace GRINS
{

  template<class Mu>
  SpalartAllmaras<Mu>::SpalartAllmaras(const std::string& physics_name, const GetPot& input )
    : TurbulenceModelsBase<Mu>(physics_name, input), // Define class variables
      _flow_vars(input,incompressible_navier_stokes),
      _turbulence_vars(input, spalart_allmaras),      
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
    this->distance_function.reset(new DistanceFunction(system->get_equation_systems(), dynamic_cast<libMesh::UnstructuredMesh&>(system->get_mesh()) ));

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
    context.get_element_fe(_turbulence_vars.nu_var())->get_JxW();
    context.get_element_fe(_turbulence_vars.nu_var())->get_phi();
    context.get_element_fe(_turbulence_vars.nu_var())->get_dphi();
    context.get_element_fe(_turbulence_vars.nu_var())->get_xyz();

    context.get_element_fe(_turbulence_vars.nu_var())->get_phi();
    context.get_element_fe(_turbulence_vars.nu_var())->get_xyz();

    context.get_side_fe(_turbulence_vars.nu_var())->get_JxW();
    context.get_side_fe(_turbulence_vars.nu_var())->get_phi();
    context.get_side_fe(_turbulence_vars.nu_var())->get_dphi();
    context.get_side_fe(_turbulence_vars.nu_var())->get_xyz();

    return;
  }

  template<class Mu>
  void SpalartAllmaras<Mu>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    const unsigned int dim = system->get_mesh().mesh_dimension();

    // Tell the system to march velocity forward in time, but
    // leave p as a constraint only
    system->time_evolving(this->_turbulence_vars.nu_var());
    
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
    libMesh::Elem &elem_pointer = context.get_elem();

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
    libMesh::AutoPtr< libMesh::DenseVector<libMesh::Real> > distance_qp;

    // Fill the vector of distances to quadrature points
    distance_qp = this->distance_function->interpolate(&elem_pointer, context.get_element_qrule().get_points());

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
        
        //const libMesh::Number x = u_qpoint[qp](0);
	//const libMesh::Number y = u_qpoint[qp](1);
	//const libMesh::Number z = (this->_dim==3)?u_qpoint[qp](2):0;

        libMesh::Real jac = JxW[qp];
	
	// The physical viscosity
	libMesh::Real _mu_qp = this->_mu(context, qp);

	// The vorticity value
	libMesh::Real _vorticity_value_qp = this->_vorticity(context, qp);
   
	// The flow velocity
	libMesh::Number u,v;
	u = context.interior_value(this->_flow_vars.u_var(), qp);
	v = context.interior_value(this->_flow_vars.v_var(), qp);
	
	libMesh::NumberVectorValue U(u,v);
	if (this->_dim == 3)
	  U(2) = context.interior_value(this->_flow_vars.w_var(), qp);
	
	//The source term
	libMesh::Real _S_tilde = this->_source_fn(nu, _mu_qp, (*distance_qp)(qp), _vorticity_value_qp);
	
	// The wall destruction term
	libMesh::Real _fw = this->_destruction_fn(nu, (*distance_qp)(qp), _S_tilde);
	
	
        // First, an i-loop over the viscosity degrees of freedom.        
        for (unsigned int i=0; i != n_nu_dofs; i++)
          {	    
	    // TODO: intialize constants cb1, cb2, cw1, sigma, and functions source_fn(nu), destruction_fn(nu), and resolve issue of grad(nu + nu_tilde)
            Fnu(i) += jac *
              ( this->_rho*(U*grad_nu)*nu_phi[i][qp]  // convection term (assumes incompressibility)
	       -this->_cb1*_S_tilde*nu*nu_phi[i][qp]  // source term
	       + (1./this->_sigma)*(-(_mu_qp+nu)*grad_nu*nu_gradphi[i][qp] - grad_nu*grad_nu*nu_phi[i][qp] + this->_cb2*grad_nu*grad_nu*nu_phi[i][qp])  // diffusion term 
               + this->_cw1*_fw*pow(nu/(*distance_qp)(qp), 2.)*nu_phi[i][qp]); // destruction term      
                    
	      // Compute the jacobian if not using numerical jacobians  
	      if (compute_jacobian)
		{
		  for (unsigned int j=0; j != n_nu_dofs; j++)
		    {
		      // TODO: precompute some terms like:
		      //   (Uvec*u_gradphi[j][qp]),
		      //   u_phi[i][qp]*u_phi[j][qp],
		      //   (u_gradphi[i][qp]*u_gradphi[j][qp])
		      
		      // Knunu(i,j) += jac *
		      // 	( 0.0     // source term
		      // 	 // diffusion term
		      // 	       ); // destruction term
                                                              
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
  libMesh::Real SpalartAllmaras<Mu>::_vorticity(AssemblyContext& context, unsigned int qp)
  {
    libMesh::Gradient grad_u, grad_v;
    grad_u = context.interior_gradient(this->_flow_vars.u_var(), qp);
    grad_v = context.interior_gradient(this->_flow_vars.v_var(), qp);

    libMesh::Real _vorticity_value;
    _vorticity_value = fabs(grad_v(0) - grad_u(1));

    if(this->_dim == 3)
      {
	libMesh::Gradient grad_w;
	grad_w = context.interior_gradient(this->_flow_vars.w_var(), qp);
	
	libMesh::Real _vorticity_component_0 = grad_w(1) - grad_v(2);
	libMesh::Real _vorticity_component_1 = grad_u(2) - grad_v(0);

	_vorticity_value += pow(pow(_vorticity_component_0, 2.0) + pow(_vorticity_component_1, 2.0) + pow(_vorticity_value, 2.0), 0.5);
      }
	return _vorticity_value;
  }
  
  template<class Mu>
    libMesh::Real SpalartAllmaras<Mu>::_source_fn(libMesh::Number nu, libMesh::Real mu, libMesh::Real wall_distance, libMesh::Real _vorticity_value)
  {
    // Step 1
    libMesh::Real _kai = nu/mu;
    
    // Step 2
    libMesh::Real _fv1 = pow(_kai, 3.0)/(pow(_kai, 3.0) + pow(this->_cv1, 3.0));

    // Step 3
    libMesh::Real _fv2 = 1 - (_kai/(1 + _kai*_fv1));

    // Step 4
    libMesh::Real _S_bar = nu/(pow(_kappa, 2.0) * pow(wall_distance, 2.0)) ;

    // Step 5, the absolute value of the vorticity
    libMesh::Real _S = _vorticity_value;

    // Step 6
    libMesh::Real _S_tilde = 0.0;
    if(_S_bar >= -this->_cv2*_S)
      {
	_S_tilde = _S + _S_bar;
      }
    else
      {
	_S_tilde = _S + (_S*(pow(this->_cv2,2.0)*_S + this->_cv3*_S_bar))/((this->_cv3 - (2*this->_cv2))*_S - _S_bar);
      }
    
    return _S_tilde;
  }

  template<class Mu>
  libMesh::Real SpalartAllmaras<Mu>::_destruction_fn(libMesh::Number nu, libMesh::Real wall_distance, libMesh::Real _S_tilde)
  {
    // Step 1
    libMesh::Real _r = 0.0;
    
    if(nu/(_S_tilde*pow(this->_kappa,2.0)*pow(wall_distance,2.0)) < this->_r_lin)
      {
	_r = nu/(_S_tilde*pow(this->_kappa,2.0)*pow(wall_distance,2.0));
      }
    else
      {
	_r = this->_r_lin;
      }

    // Step 2
    libMesh::Real _g = _r + this->_c_w2*(pow(_r,6.0) - _r);

    // Step 3
    libMesh::Real _fw = _g*pow((1 + pow(this->_c_w3,6.0))/(pow(_g,6.0) + pow(this->_c_w3,6.0)), 1.0/6.0);
    
    return _fw;
  }
  

  
} // namespace GRINS

// Instantiate
INSTANTIATE_TURBULENCE_MODELS_SUBCLASS(SpalartAllmaras);
