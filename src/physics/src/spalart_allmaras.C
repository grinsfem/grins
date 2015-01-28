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
#include "grins/spalart_allmaras_bc_handling.h"
#include "grins/turbulence_models_macro.h"

#include "grins/constant_viscosity.h"
#include "grins/parsed_viscosity.h"
#include "grins/spalart_allmaras_viscosity.h"

// libMesh
#include "libmesh/quadrature.h"
#include "libmesh/elem.h"

namespace GRINS
{

  template<class Mu>
  SpalartAllmaras<Mu>::SpalartAllmaras(const std::string& physics_name, const GetPot& input )
    : TurbulenceModelsBase<Mu>(physics_name, input), // Define class variables      
      _flow_vars(input,incompressible_navier_stokes),
      _turbulence_vars(input, spalart_allmaras),      
      _spalart_allmaras_helper(input),
      _no_of_walls(input("Physics/"+spalart_allmaras+"/no_of_walls", 0)),
      boundary_mesh(libMesh::NULL, _dim)
  {    
    // Loop over the _no_of_walls and fill the wall_ids set
    for(unsigned int i = 0; i != _no_of_walls; i++)
      _wall_ids.insert(input("Physics/"+spalart_allmaras+"/wall_ids", 0.0, i ));
    
    // This is deleted in the base class
    this->_bc_handler = new SpalartAllmarasBCHandling( physics_name, input );
    
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
    // Init base class.
    TurbulenceModelsBase<Mu>::init_variables(system);

    // Use the _wall_ids set to build the boundary mesh object
    (system->get_mesh()).boundary_info->sync(_wall_ids, boundary_mesh);        

    //this->distance_function.reset(new DistanceFunction(system->get_equation_systems(), dynamic_cast<libMesh::UnstructuredMesh&>(system->get_mesh()) ));
    this->distance_function.reset(new DistanceFunction(system->get_equation_systems(), boundary_mesh));
                 
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
        
        const libMesh::Number x = nu_qpoint[qp](0);
	const libMesh::Number y = nu_qpoint[qp](1);
	//const libMesh::Number z = (this->_dim==3)?nu_qpoint[qp](2):0;

        libMesh::Real jac = JxW[qp];
	
	// The physical viscosity	
	libMesh::Real _mu_qp = this->_mu(context, qp);

	// The vorticity value
	libMesh::Real _vorticity_value_qp = this->_spalart_allmaras_helper._vorticity(context, qp);
   
	// The flow velocity
	libMesh::Number u,v;
	u = context.interior_value(this->_flow_vars.u_var(), qp);
	v = context.interior_value(this->_flow_vars.v_var(), qp);
	
	libMesh::NumberVectorValue U(u,v);
	if (this->_dim == 3)
	  U(2) = context.interior_value(this->_flow_vars.w_var(), qp);
	
	// The calculated distance
	std::cout<<"Distance to wall from point("<<x<<","<<y<<") is: "<< ( (*distance_qp)(qp) ) <<std::endl;

	//The source term
	libMesh::Real _S_tilde = this->_spalart_allmaras_helper._source_fn(nu, _mu_qp, (*distance_qp)(qp), _vorticity_value_qp);
	libMesh::Real _source_term = ((*distance_qp)(qp)==0.0)?1.0:this->_spalart_allmaras_helper._cb1*_S_tilde*nu;

	// The wall destruction term
	libMesh::Real _fw = this->_spalart_allmaras_helper._destruction_fn(nu, (*distance_qp)(qp), _S_tilde);
	libMesh::Real _destruction_term = ((*distance_qp)(qp)==0.0)?1.0:this->_spalart_allmaras_helper._cw1*_fw*pow(nu/(*distance_qp)(qp), 2.);
	
        // First, an i-loop over the viscosity degrees of freedom.        
        for (unsigned int i=0; i != n_nu_dofs; i++)
          {	    
	    // TODO: intialize constants cb1, cb2, cw1, sigma, and functions source_fn(nu), destruction_fn(nu), and resolve issue of grad(nu + nu_tilde)
            Fnu(i) += jac *
              ( this->_rho*(U*grad_nu)*nu_phi[i][qp]  // convection term (assumes incompressibility)
	    -_source_term*nu_phi[i][qp] // source term
	      + (1./this->_spalart_allmaras_helper._sigma)*(-(_mu_qp+nu)*grad_nu*nu_gradphi[i][qp] - grad_nu*grad_nu*nu_phi[i][qp] + this->_spalart_allmaras_helper._cb2*grad_nu*grad_nu*nu_phi[i][qp]) // diffusion term 
		+ _destruction_term*nu_phi[i][qp]); // destruction term
	    
	    //Fnu(i) += jac * (grad_nu*nu_gradphi[i][qp]);
                    
	      // Compute the jacobian if not using numerical jacobians  
	      if (compute_jacobian)
		{
		  libmesh_not_implemented();		     
		} // end - if (compute_jacobian)

          } // end of the outer dof (i) loop
      } // end of the quadrature point (qp) loop

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("SpalartAllmaras::element_time_derivative");
#endif

    return;
  }    
  
} // namespace GRINS

// Instantiate
INSTANTIATE_TURBULENCE_MODELS_SUBCLASS(SpalartAllmaras);
