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

// This class
#include "reacting_low_mach_navier_stokes.h"

// GRINS
#include "cached_quantities_enum.h"

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
	    /*! \todo Need to figure out something smarter for controling species
	              that go slightly negative. */
	    const Real value = std::max( c.interior_value(this->_species_vars[s],qp), 0.0 );
	    libmesh_assert_greater_equal(value,0.0);
	    Y.push_back(value);
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
	const Real T = c.interior_value(this->_T_var, qp);
	libmesh_assert_greater(T, 0.0);


	const Real p0 = this->get_p0_steady(c,qp);
	libmesh_assert_greater(p0, 0.0);

	ReactingFlowCache cache( T, p0, Y );

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
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

    const GRINS::BoundaryID boundary_id =
      system->get_mesh().boundary_info->boundary_id(c.elem, c.side);

    libmesh_assert (boundary_id != libMesh::BoundaryInfo::invalid_id);

    this->_bc_handler->apply_neumann_bcs( c, this->_species_vars[3], request_jacobian, boundary_id );

    this->_bc_handler->apply_neumann_bcs( c, this->_species_vars[1], request_jacobian, boundary_id );
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
    
    libMesh::DenseSubVector<Number>& Fp = *c.elem_subresiduals[this->_p_var]; // R_{p}
    
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
    
    const libMesh::Number term1 = -U*(mass_term + grad_T/T);

    const libMesh::Number termf = (term1 + divU)*JxW[qp];
      
    for (unsigned int i=0; i != n_p_dofs; i++)
      {
	Fp(i) += termf*p_phi[i][qp];
	libmesh_assert( !libmesh_isnan(Fp(i)) );
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

	const Real term1 = -rho*(U*grad_w[s]) + omega_dot[s];
	const libMesh::Gradient term2 = -rho*D[s]*grad_w[s];

	for (unsigned int i=0; i != n_s_dofs; i++)
	  {
	    /*! \todo Need to add SCEBD term. */
	    Fs(i) += ( term1*s_phi[i][qp] + term2*s_grad_phi[i][qp] )*JxW[qp];

	    libmesh_assert( !libmesh_isnan(Fs(i)) );
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

	libmesh_assert( !libmesh_isnan(Fu(i)) );

	Fv(i) += ( -rho*U*grad_v*u_phi[i][qp]                 // convection term
		   + p*u_gradphi[i][qp](1)                           // pressure term
		   - mu*(u_gradphi[i][qp]*grad_v + u_gradphi[i][qp]*grad_vT
			 - 2.0/3.0*divU*u_gradphi[i][qp](1) )    // diffusion term
		   + rho*this->_g(1)*u_phi[i][qp]                 // hydrostatic term
		   )*JxW[qp];

	libmesh_assert( !libmesh_isnan(Fv(i)) );

	if (this->_dim == 3)
	  {
	    Fw(i) += ( -rho*U*grad_w*u_phi[i][qp]                 // convection term
		       + p*u_gradphi[i][qp](2)                           // pressure term
		       - mu*(u_gradphi[i][qp]*grad_w + u_gradphi[i][qp]*grad_wT
			     - 2.0/3.0*divU*u_gradphi[i][qp](2) )    // diffusion term
		       + rho*this->_g(2)*u_phi[i][qp]                 // hydrostatic term
		       )*JxW[qp];

	    libmesh_assert( !libmesh_isnan(Fw(i)) );
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

    libmesh_assert( !libmesh_isnan(chem_term) );

    for (unsigned int i=0; i != n_T_dofs; i++)
      {
	FT(i) += ( ( -rho*cp*U*grad_T + chem_term )*T_phi[i][qp] // convection term + chemistry term
		     - k*grad_T*T_gradphi[i][qp]   /* diffusion term */   )*JxW[qp];

	libmesh_assert( !libmesh_isnan(FT(i)) );
      }

    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::compute_cache( const libMesh::FEMContext& context, 
							    CachedValues& cache )
  {
    libmesh_not_implemented();
    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::compute_cache( const libMesh::FEMContext& context, 
							    const std::vector<libMesh::Point>& points,
							    CachedValues& cache )
  {
    if( cache.is_active(Cache::MIXTURE_DENSITY) )
      {
	std::vector<Real> rho_values;
	rho_values.reserve( points.size() );

	std::vector<Real> mass_fracs( this->_n_species );
	
	for( std::vector<libMesh::Point>::const_iterator point = points.begin();
	     point != points.end(); point++ )
	  {
	    Real T = this->T(*point,context);
	    Real p0 = this->get_p0_steady(context,*point);
	    this->mass_fractions( *point, context, mass_fracs );

	    rho_values.push_back(this->rho( T, p0, mass_fracs) );
	  }

	cache.set_values( Cache::MIXTURE_DENSITY, rho_values );
      }

    if( cache.is_active(Cache::MOLE_FRACTIONS) )
      {
	std::vector<std::vector<Real> > mole_fractions;
	mole_fractions.resize( points.size() );

	for( unsigned int i = 0; i < points.size(); i++ )
	  {
	    mole_fractions[i].resize( this->_n_species );
	  }

	std::vector<Real> mass_fracs( this->_n_species );

	for( unsigned int p = 0; p < points.size(); p++ )
	  {
	    this->mass_fractions( points[p], context, mass_fracs );

	    Real M = this->_gas_mixture.M(mass_fracs);

	    this->_gas_mixture.X(M,mass_fracs,mole_fractions[p]);
	  }

	cache.set_vector_values(Cache::MOLE_FRACTIONS, mole_fractions );
      }

    return;
  }

  // Instantiate
#ifdef GRINS_HAVE_CANTERA
  template class ReactingLowMachNavierStokes< IdealGasMixture<CanteraThermodynamics,CanteraTransport,CanteraKinetics> >;
  template class ReactingLowMachNavierStokes< IdealGasMixture<CanteraThermodynamics,ConstantTransport,CanteraKinetics> >;
#endif

} // namespace GRINS
