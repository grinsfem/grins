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
#include "grins/spalart_allmaras_spgsm_stab.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/constant_viscosity.h"
#include "grins/parsed_viscosity.h"
#include "grins/spalart_allmaras_viscosity.h"
#include "grins/turbulence_models_macro.h"

//libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{

  template<class Mu>
  SpalartAllmarasSPGSMStabilization<Mu>::SpalartAllmarasSPGSMStabilization( const std::string& physics_name, 
                                                                                              const GetPot& input )
    : SpalartAllmarasStabilizationBase<Mu>(physics_name,input)
  {
    this->read_input_options(input);

    return;
  }

  template<class Mu>
  SpalartAllmarasSPGSMStabilization<Mu>::~SpalartAllmarasSPGSMStabilization()
  {
    return;
  }
  
  template<class Mu>
  void SpalartAllmarasSPGSMStabilization<Mu>::element_time_derivative( bool compute_jacobian,
                                                                              AssemblyContext& context,
                                                                              CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("SpalartAllmarasSPGSMStabilization::element_time_derivative");
#endif

    // Get a pointer to the current element, we need this for computing the distance to wall for the
    // quadrature points
    libMesh::Elem &elem_pointer = context.get_elem();

    // The number of local degrees of freedom in each variable.
    const unsigned int n_nu_dofs = context.get_dof_indices(this->_turbulence_vars.nu_var()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_turbulence_vars.u_var())->get_JxW();

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
    
    libMesh::DenseSubMatrix<libMesh::Number> &Knunu = context.get_elem_jacobian(this->_turbulence_vars.nu_var(), this->_turbulence_vars.nu_var()); // R_{nu},{nu}
    
    libMesh::DenseSubVector<libMesh::Number> &Fnu = context.get_elem_residual(this->_turbulence_vars.nu_var()); // R_{nu}
    
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

	// Stabilization terms

	libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );
			
        libMesh::Real tau_spalart = this->_stab_helper.compute_tau_spalart( context, qp, g, G, this->_rho, U, _mu_qp, this->_is_steady );
        
        libMesh::Number RM_spalart = this->_stab_helper.compute_res_spalart_steady( context, qp, this->_rho, _mu_qp, (*distance_qp)(qp) );
        
        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fnu(i) += jac*( tau_spalart*RM_spalart*nu_phi[i][qp]  );            
          }

        if( compute_jacobian )
          {
            libmesh_not_implemented();
          }

      }

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("SpalartAllmarasSPGSMStabilization::element_time_derivative");
#endif

    return;
  }

  template<class Mu>
  void SpalartAllmarasSPGSMStabilization<Mu>::element_constraint( bool compute_jacobian,
                                                                         AssemblyContext& context,
                                                                         CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("SpalartAllmarasSPGSMStabilization::element_constraint");
#endif

    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.get_dof_indices(this->_flow_vars.p_var()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u_var())->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
      context.get_element_fe(this->_flow_vars.p_var())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(this->_flow_vars.p_var()); // R_{p}

    libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u_var());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );
      
        libMesh::RealGradient U( context.interior_value( this->_flow_vars.u_var(), qp ),
                                 context.interior_value( this->_flow_vars.v_var(), qp ) );
        if( this->_dim == 3 )
          {
            U(2) = context.interior_value( this->_flow_vars.w_var(), qp );
          }

	// Compute the viscosity at this qp
	libMesh::Real _mu_qp = this->_mu(context, qp);

        libMesh::Real tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, this->_rho, U, _mu_qp, this->_is_steady );

        libMesh::RealGradient RM_s = this->_stab_helper.compute_res_momentum_steady( context, qp, this->_rho, _mu_qp );

        for (unsigned int i=0; i != n_p_dofs; i++)
          {
            Fp(i) += tau_M*RM_s*p_dphi[i][qp]*JxW[qp];
          }

        if( compute_jacobian )
          {
            libmesh_not_implemented();
          }

      }

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("SpalartAllmarasSPGSMStabilization::element_constraint");
#endif

    return;
  }

  template<class Mu>
  void SpalartAllmarasSPGSMStabilization<Mu>::mass_residual( bool compute_jacobian,
                                                                    AssemblyContext& context,
                                                                    CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("SpalartAllmarasSPGSMStabilization::mass_residual");
#endif

    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.get_dof_indices(this->_flow_vars.p_var()).size();
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u_var()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u_var())->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
      context.get_element_fe(this->_flow_vars.p_var())->get_dphi();

    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_flow_vars.u_var())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u_var()); // R_{p}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v_var()); // R_{p}
    libMesh::DenseSubVector<libMesh::Number> *Fw = NULL;
    if(this->_dim == 3)
      Fw = &context.get_elem_residual(this->_flow_vars.w_var()); // R_{w}

    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(this->_flow_vars.p_var()); // R_{p}

    libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u_var());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::RealGradient U( context.fixed_interior_value( this->_flow_vars.u_var(), qp ),
                                 context.fixed_interior_value( this->_flow_vars.v_var(), qp ) );
	// Compute the viscosity at this qp
	libMesh::Real _mu_qp = this->_mu(context, qp);

        if( this->_dim == 3 )
          {
            U(2) = context.fixed_interior_value( this->_flow_vars.w_var(), qp );
          }
      
        libMesh::Real tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, this->_rho, U, _mu_qp, false );

        libMesh::RealGradient RM_t = this->_stab_helper.compute_res_momentum_transient( context, qp, this->_rho );

        for (unsigned int i=0; i != n_p_dofs; i++)
          {
            Fp(i) -= tau_M*RM_t*p_dphi[i][qp]*JxW[qp];
          }

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) -= tau_M*RM_t(0)*this->_rho*U*u_gradphi[i][qp]*JxW[qp];

            Fv(i) -= tau_M*RM_t(1)*this->_rho*U*u_gradphi[i][qp]*JxW[qp];

            if( this->_dim == 3 )
              {
                (*Fw)(i) -= tau_M*RM_t(2)*this->_rho*U*u_gradphi[i][qp]*JxW[qp];
              }
          }

        if( compute_jacobian )
          {
            libmesh_not_implemented();
          }

      }

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("SpalartAllmarasSPGSMStabilization::mass_residual");
#endif

    return;
  }

} // end namespace GRINS

// Instantiate
INSTANTIATE_TURBULENCE_MODELS_SUBCLASS(SpalartAllmarasSPGSMStabilization);
