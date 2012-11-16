//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "reacting_low_mach_navier_stokes.h"

namespace GRINS
{
  template<class Mixture>
  ReactingLowMachNavierStokes<Mixture>::ReactingLowMachNavierStokes(const PhysicsName& physics_name, const GetPot& input)
    : ReactingLowMachNavierStokesBase<Mixture>(physics_name,input),
      _p_pinning(input,physics_name)
  {
    this->read_input_options(input);

    // This is deleted in the base class
    this->_bc_handler = new ReactingLowMachNavierStokesBCHandling( physics_name, input );

    return;
  }

  template<class Mixture>
  ReactingLowMachNavierStokes<Mixture>::~ReactingLowMachNavierStokes()
  {
    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::read_input_options( const GetPot& input )
  {
    // Other quantities read in base class

    // Read pressure pinning information
    this->_pin_pressure = input("Physics/"+reacting_low_mach_navier_stokes+"/pin_pressure", false );
  
    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::init_context( libMesh::DiffContext &context )
  {
    // First call base class
    GRINS::ReactingLowMachNavierStokesBase<Mixture>::init_context(context);

    libMesh::FEMContext &c = libmesh_cast_ref<libMesh::FEMContext&>(context);

    // We also need the side shape functions, etc.
    c.side_fe_var[this->_u_var]->get_JxW();
    c.side_fe_var[this->_u_var]->get_phi();
    c.side_fe_var[this->_u_var]->get_dphi();
    c.side_fe_var[this->_u_var]->get_xyz();

    c.side_fe_var[this->_T_var]->get_JxW();
    c.side_fe_var[this->_T_var]->get_phi();
    c.side_fe_var[this->_T_var]->get_dphi();
    c.side_fe_var[this->_T_var]->get_xyz();

    return;
  }

  template<class Mixture>
  bool ReactingLowMachNavierStokes<Mixture>::element_time_derivative( bool request_jacobian,
								      libMesh::DiffContext& context,
								      libMesh::FEMSystem* system )
  {
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    
    unsigned int n_qpoints = c.element_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	// Compute species mass fractions at quadrature points
	std::vector<Real> Y;
	Y.reserve(this->_n_species);
	for( unsigned int s = 0; s < this->_n_species; s++ )
	  {
	    Y.push_back( c.interior_value(this->_species_vars[s],qp) );
	  }

	// Build up cache at element interior quadrature point
	/*! \todo Ought to rethink constructing this at every quadrature point. Perhaps
	          add a "reset" method (or something similar) that just clears everything
	          so we don't have to keep deallocating/reallocating. Probably will require
	          adding libmesh_asserts in the themro/transport/chemistry calculations since
	          we assume at least T, P, and Y are already present in the cache.
	          Actually, what we want to do is have the cache built once per element
	          (thus caching at every quadrature point), so this can be reused for multiple
	          Physics classes. */
	ReactingFlowCache cache( c.interior_value(this->_T_var, qp), 
				 this->get_p0_steady(c,qp), Y );

	this->build_reacting_flow_cache(c, cache, qp);

	this->assemble_mass_time_deriv(c, cache, qp);
	this->assemble_species_time_deriv(c, cache, qp);
	this->assemble_momentum_time_deriv(c, cache, qp);
	this->assemble_energy_time_deriv(c, cache, qp);
      }

    // Pin p = p_value at p_point
    if( this->_pin_pressure )
      {
	this->_p_pinning.pin_value( context, request_jacobian, this->_p_var );
      }

    return request_jacobian;
  }


  template<class Mixture>
  bool ReactingLowMachNavierStokes<Mixture>::element_constraint( bool request_jacobian,
								 libMesh::DiffContext& context,
								 libMesh::FEMSystem* system )
  {
#ifdef USE_GRVY_TIMERS
    //this->_timer->BeginTimer("LowMachNavierStokes::element_constraint");
#endif

    //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

#ifdef USE_GRVY_TIMERS
    //this->_timer->EndTimer("LowMachNavierStokes::element_constraint");
#endif

    return request_jacobian;
  }

  template<class Mixture>
  bool ReactingLowMachNavierStokes<Mixture>::side_time_derivative( bool request_jacobian,
								   libMesh::DiffContext& context,
								   libMesh::FEMSystem* system )
  {
    /*! \todo Need to implement thermodynamic pressure calcuation for cases where it's needed. */
    /*
    if( this->_enable_thermo_press_calc )
      {
#ifdef USE_GRVY_TIMERS
	this->_timer->BeginTimer("LowMachNavierStokes::side_time_derivative");
#endif
	FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

	this->assemble_thermo_press_side_time_deriv( request_jacobian, c, system );

#ifdef USE_GRVY_TIMERS
	this->_timer->EndTimer("LowMachNavierStokes::side_time_derivative");
#endif
      }
    */
    return request_jacobian;
  }

  template<class Mixture>
  bool ReactingLowMachNavierStokes<Mixture>::side_constraint( bool request_jacobian,
								libMesh::DiffContext& context,
								libMesh::FEMSystem* system )
  {
#ifdef USE_GRVY_TIMERS
    //this->_timer->BeginTimer("ReactingLowMachNavierStokes::side_constraint");
#endif

    //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

#ifdef USE_GRVY_TIMERS
    //this->_timer->EndTimer("ReactingLowMachNavierStokes::side_constraint");
#endif

    return request_jacobian;
  }

  template<class Mixture>
  bool ReactingLowMachNavierStokes<Mixture>::mass_residual( bool request_jacobian,
							    libMesh::DiffContext& context,
							    libMesh::FEMSystem* system )
  {
    libmesh_not_implemented();
    /*
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    */
    return request_jacobian;
  }


  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::assemble_mass_time_deriv(libMesh::FEMContext& c, 
								      const ReactingFlowCache& cache, 
								      unsigned int qp)
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = c.dof_indices_var[this->_p_var].size();
    
    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      c.element_fe_var[this->_u_var]->get_JxW();
    
    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& p_phi =
      c.element_fe_var[this->_p_var]->get_phi();
    
    libMesh::DenseSubVector<Number> &Fp = *c.elem_subresiduals[this->_p_var]; // R_{p}
    
    libMesh::Number T = cache.T();
    const libMesh::NumberVectorValue& U = cache.U();
    libMesh::Number divU = cache.divU();
    const libMesh::Gradient& grad_T = cache.grad_T();
    Real M = cache.M();
    const std::vector<libMesh::Gradient>& grad_w = cache.mass_fractions_grad();

    libmesh_assert_equal_to( grad_w.size(), this->_n_species );
    
    libMesh::Gradient mass_term(0.0,0.0,0.0);
    for(unsigned int s=0; s < this->_n_species; s++ )
      {
	mass_term += grad_w[s]/this->_gas_mixture.M(s);
      }
    mass_term *= M;
    
    for (unsigned int i=0; i != n_p_dofs; i++)
      {
	Fp(i) += (-U*(mass_term + grad_T/T) + divU)*p_phi[i][qp]*JxW[qp];
      }

    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::assemble_species_time_deriv(libMesh::FEMContext& c, 
									 const ReactingFlowCache& cache, 
									 unsigned int qp)
  {
    // Convenience
    const VariableIndex s0_var = this->_species_vars[0];
    
    /* The number of local degrees of freedom in each species variable.
       We assume the same number of dofs for each species */
    const unsigned int n_s_dofs = c.dof_indices_var[s0_var].size();
    
    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW = c.element_fe_var[s0_var]->get_JxW();
    
    // The species shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& s_phi = c.element_fe_var[s0_var]->get_phi();

    // The species shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::Gradient> >& s_grad_phi = c.element_fe_var[s0_var]->get_dphi();

    libMesh::Number rho = cache.rho();
    const libMesh::NumberVectorValue& U = cache.U();
    const std::vector<libMesh::Gradient>& grad_w = cache.mass_fractions_grad();
    const std::vector<Real>& D = cache.D();
    const std::vector<Real>& omega_dot = cache.omega_dot();

    for(unsigned int s=0; s < this->_n_species; s++ )
      {
	libMesh::DenseSubVector<Number> &Fs = *c.elem_subresiduals[this->_species_vars[s]]; // R_{s}

	for (unsigned int i=0; i != n_s_dofs; i++)
	  {
	    /*! \todo Need to add SCEBD term. */
	    Fs(i) += ( ( rho*(U*grad_w[s]) - omega_dot[s] )*s_phi[i][qp] 
		       + rho*D[s]*(grad_w[s]*s_grad_phi[i][qp])        )*JxW[qp];
	  }
      }

    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::assemble_momentum_time_deriv(libMesh::FEMContext& c, 
									  const ReactingFlowCache& cache, 
									  unsigned int qp)
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = c.dof_indices_var[this->_u_var].size();

    // Check number of dofs is same for _u_var, v_var and w_var.
    libmesh_assert (n_u_dofs == c.dof_indices_var[this->_v_var].size());
    if (this->_dim == 3)
      libmesh_assert (n_u_dofs == c.dof_indices_var[this->_w_var].size());

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      c.element_fe_var[this->_u_var]->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& u_phi =
      c.element_fe_var[this->_u_var]->get_phi();

    // The velocity shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      c.element_fe_var[this->_u_var]->get_dphi();

    libMesh::DenseSubVector<Number> &Fu = *c.elem_subresiduals[this->_u_var]; // R_{u}
    libMesh::DenseSubVector<Number> &Fv = *c.elem_subresiduals[this->_v_var]; // R_{v}
    libMesh::DenseSubVector<Number> &Fw = *c.elem_subresiduals[this->_w_var]; // R_{w}
    
    libMesh::Number rho = cache.rho();
    const libMesh::NumberVectorValue& U = cache.U();
    libMesh::Number divU = cache.divU();
    libMesh::Number mu = cache.mu();
    libMesh::Number p = cache.p_hydro();
    const libMesh::Gradient& grad_u = cache.grad_u();
    const libMesh::Gradient& grad_v = cache.grad_v();
    const libMesh::Gradient& grad_w = cache.grad_w();

    libMesh::NumberVectorValue grad_uT( grad_u(0), grad_v(0) ); 
    libMesh::NumberVectorValue grad_vT( grad_u(1), grad_v(1) );
    libMesh::NumberVectorValue grad_wT;
    if( this->_dim == 3 )
      {
	grad_uT(2) = grad_w(0);
	grad_vT(2) = grad_w(1);
	grad_wT = libMesh::NumberVectorValue( grad_u(2), grad_v(2), grad_w(2) );
      }

    // Now a loop over the pressure degrees of freedom.  This
    // computes the contributions of the continuity equation.
    for (unsigned int i=0; i != n_u_dofs; i++)
      {
	Fu(i) += ( -rho*U*grad_u*u_phi[i][qp]                 // convection term
		   + p*u_gradphi[i][qp](0)                           // pressure term
		   - mu*(u_gradphi[i][qp]*grad_u + u_gradphi[i][qp]*grad_uT
			 - 2.0/3.0*divU*u_gradphi[i][qp](0) )    // diffusion term
		   + rho*this->_g(0)*u_phi[i][qp]                 // hydrostatic term
		   )*JxW[qp]; 

	Fv(i) += ( -rho*U*grad_v*u_phi[i][qp]                 // convection term
		   + p*u_gradphi[i][qp](1)                           // pressure term
		   - mu*(u_gradphi[i][qp]*grad_v + u_gradphi[i][qp]*grad_vT
			 - 2.0/3.0*divU*u_gradphi[i][qp](1) )    // diffusion term
		   + rho*this->_g(1)*u_phi[i][qp]                 // hydrostatic term
		   )*JxW[qp];

	if (this->_dim == 3)
	  {
	    Fw(i) += ( -rho*U*grad_w*u_phi[i][qp]                 // convection term
		       + p*u_gradphi[i][qp](2)                           // pressure term
		       - mu*(u_gradphi[i][qp]*grad_w + u_gradphi[i][qp]*grad_wT
			     - 2.0/3.0*divU*u_gradphi[i][qp](2) )    // diffusion term
		       + rho*this->_g(2)*u_phi[i][qp]                 // hydrostatic term
		       )*JxW[qp];
	  }
      }
    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::assemble_energy_time_deriv(libMesh::FEMContext& c, 
									const ReactingFlowCache& cache, 
									unsigned int qp)
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = c.dof_indices_var[this->_T_var].size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      c.element_fe_var[this->_T_var]->get_JxW();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      c.element_fe_var[this->_T_var]->get_phi();

    // The temperature shape functions gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      c.element_fe_var[this->_T_var]->get_dphi();

    libMesh::DenseSubVector<Number> &FT = *c.elem_subresiduals[this->_T_var]; // R_{T}

    libMesh::Number rho = cache.rho() ;
    const libMesh::NumberVectorValue& U = cache.U();
    libMesh::Number cp = cache.cp();
    libMesh::Number k = cache.k();
    const libMesh::Gradient& grad_T = cache.grad_T();
    const std::vector<Real>& omega_dot = cache.omega_dot();
    const std::vector<Real>& h = cache.species_enthalpy();

    Real chem_term = 0.0;
    for(unsigned int s=0; s < this->_n_species; s++ )
      {
	chem_term += h[s]*omega_dot[s];
      }

    for (unsigned int i=0; i != n_T_dofs; i++)
      {
	FT(i) += ( ( -rho*cp*U*grad_T + chem_term )*T_phi[i][qp] // convection term + chemistry term
		     - k*grad_T*T_gradphi[i][qp]   /* diffusion term */   )*JxW[qp]; 
      }

    return;
  }

  // Instantiate
#ifdef HAVE_CANTERA
  template class ReactingLowMachNavierStokes< IdealGasMixture<CanteraThermodynamics,CanteraTransport,CanteraKinetics> >;
#endif

} // namespace GRINS
