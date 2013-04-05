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
#include "grins/reacting_low_mach_navier_stokes.h"

// GRINS
#include "grins/cached_quantities_enum.h"
#include "grins/cea_thermo.h"
#include "grins/cantera_thermo.h"
#include "grins/constant_transport.h"
#include "grins/cantera_transport.h"
#include "grins/cantera_kinetics.h"
#include "grins/grins_kinetics.h"
#include "grins/reacting_low_mach_navier_stokes_bc_handling.h"

// libMesh
#include "libmesh/quadrature.h"
#include "libmesh/fem_system.h"
#include "libmesh/fem_context.h"

namespace GRINS
{
  template<class Mixture>
  ReactingLowMachNavierStokes<Mixture>::ReactingLowMachNavierStokes(const PhysicsName& physics_name, const GetPot& input)
    : ReactingLowMachNavierStokesBase<Mixture>(physics_name,input),
      _p_pinning(input,physics_name)
  {
    this->read_input_options(input);

    // This is deleted in the base class
    this->_bc_handler = new ReactingLowMachNavierStokesBCHandling( physics_name, input,
								   this->_gas_mixture.chem_mixture() );

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
  void ReactingLowMachNavierStokes<Mixture>::init_context( libMesh::FEMContext& context )
  {
    // First call base class
    GRINS::ReactingLowMachNavierStokesBase<Mixture>::init_context(context);

    // We also need the side shape functions, etc.
    context.side_fe_var[this->_u_var]->get_JxW();
    context.side_fe_var[this->_u_var]->get_phi();
    context.side_fe_var[this->_u_var]->get_dphi();
    context.side_fe_var[this->_u_var]->get_xyz();

    context.side_fe_var[this->_T_var]->get_JxW();
    context.side_fe_var[this->_T_var]->get_phi();
    context.side_fe_var[this->_T_var]->get_dphi();
    context.side_fe_var[this->_T_var]->get_xyz();

    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::element_time_derivative( bool compute_jacobian,
								      libMesh::FEMContext& context,
								      CachedValues& cache )
  {
    unsigned int n_qpoints = context.element_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	this->assemble_mass_time_deriv(context, qp, cache);
	this->assemble_species_time_deriv(context, qp, cache);
	this->assemble_momentum_time_deriv(context, qp, cache);
	this->assemble_energy_time_deriv(context, qp, cache);
      }

    // Pin p = p_value at p_point
    if( this->_pin_pressure )
      {
	this->_p_pinning.pin_value( context, compute_jacobian, this->_p_var );
      }

    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::side_time_derivative( bool compute_jacobian,
								   libMesh::FEMContext& context,
								   CachedValues& cache )
  {
    /*! \todo Need to implement thermodynamic pressure calcuation for cases where it's needed. */

    std::vector<BoundaryID> ids = context.side_boundary_ids();

    for( std::vector<BoundaryID>::const_iterator it = ids.begin();
	 it != ids.end(); it++ )
      {
	libmesh_assert (*it != libMesh::BoundaryInfo::invalid_id);

	this->_bc_handler->apply_neumann_bcs( context, cache, compute_jacobian, *it );
      }

    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::mass_residual( bool /*compute_jacobian*/,
							    libMesh::FEMContext& /*context*/,
							    CachedValues& /*cache*/ )
  {
    libmesh_not_implemented();
    /*
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    */
    return;
  }


  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::assemble_mass_time_deriv( libMesh::FEMContext& context, 
								       unsigned int qp,
								       const CachedValues& cache )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.dof_indices_var[this->_p_var].size();
    
    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real>& JxW =
      context.element_fe_var[this->_u_var]->get_JxW();
    
    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& p_phi =
      context.element_fe_var[this->_p_var]->get_phi();
    
    const std::vector<libMesh::Point>& u_qpoint = 
      context.element_fe_var[this->_u_var]->get_xyz();

    libMesh::DenseSubVector<libMesh::Number>& Fp = *context.elem_subresiduals[this->_p_var]; // R_{p}

    libMesh::Number u, v, w, T;
    u = cache.get_cached_values(Cache::X_VELOCITY)[qp];
    v = cache.get_cached_values(Cache::Y_VELOCITY)[qp];
    if (this->_dim == 3)
      w = cache.get_cached_values(Cache::Z_VELOCITY)[qp];
    
    T = cache.get_cached_values(Cache::TEMPERATURE)[qp];

    libMesh::NumberVectorValue U(u,v);
    if (this->_dim == 3)
      U(2) = w;

    libMesh::Gradient grad_u = cache.get_cached_gradient_values(Cache::X_VELOCITY_GRAD)[qp];
    libMesh::Gradient grad_v = cache.get_cached_gradient_values(Cache::Y_VELOCITY_GRAD)[qp];

    libMesh::Gradient grad_w;
    if (this->_dim == 3)
      grad_w = cache.get_cached_gradient_values(Cache::Z_VELOCITY_GRAD)[qp];

    libMesh::Number divU = grad_u(0) + grad_v(1);
    if (this->_dim == 3)
      divU += grad_w(2);

    const libMesh::Number r = u_qpoint[qp](0);

    libMesh::Real jac = JxW[qp];

    if( this->_bc_handler->is_axisymmetric() )
      {
	divU += U(0)/r;
	jac *= r;
      }

    libMesh::Gradient grad_T = 
      cache.get_cached_gradient_values(Cache::TEMPERATURE_GRAD)[qp];

    libMesh::Real M = cache.get_cached_values(Cache::MOLAR_MASS)[qp];

    std::vector<libMesh::Gradient> grad_ws = cache.get_cached_vector_gradient_values(Cache::MASS_FRACTIONS_GRAD)[qp];
    libmesh_assert_equal_to( grad_ws.size(), this->_n_species );
    
    libMesh::Gradient mass_term(0.0,0.0,0.0);
    for(unsigned int s=0; s < this->_n_species; s++ )
      {
	mass_term += grad_ws[s]/this->_gas_mixture.M(s);
      }
    mass_term *= M;
    
    const libMesh::Number term1 = -U*(mass_term + grad_T/T);

    const libMesh::Number termf = (term1 + divU)*jac;
      
    for (unsigned int i=0; i != n_p_dofs; i++)
      {
	Fp(i) += termf*p_phi[i][qp];
	libmesh_assert( !libmesh_isnan(Fp(i)) );
      }

    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::assemble_species_time_deriv(libMesh::FEMContext& context, 
									 unsigned int qp,
									 const CachedValues& cache)
  {
    // Convenience
    const VariableIndex s0_var = this->_species_vars[0];
    
    /* The number of local degrees of freedom in each species variable.
       We assume the same number of dofs for each species */
    const unsigned int n_s_dofs = context.dof_indices_var[s0_var].size();
    
    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW = context.element_fe_var[s0_var]->get_JxW();
    
    // The species shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& s_phi = context.element_fe_var[s0_var]->get_phi();

    // The species shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::Gradient> >& s_grad_phi = context.element_fe_var[s0_var]->get_dphi();

    const std::vector<libMesh::Point>& s_qpoint = 
      context.element_fe_var[this->_species_vars[0]]->get_xyz();

    libMesh::Number rho = cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp];

    libMesh::Number u, v, w;
    u = cache.get_cached_values(Cache::X_VELOCITY)[qp];
    v = cache.get_cached_values(Cache::Y_VELOCITY)[qp];
    if (this->_dim == 3)
      w = cache.get_cached_values(Cache::Z_VELOCITY)[qp];

    libMesh::NumberVectorValue U(u,v);
    if (this->_dim == 3)
      U(2) = w;

    std::vector<libMesh::Gradient> grad_w = 
      cache.get_cached_vector_gradient_values(Cache::MASS_FRACTIONS_GRAD)[qp];
    libmesh_assert_equal_to( grad_w.size(), this->_n_species );

    const std::vector<libMesh::Real>& D = 
      cache.get_cached_vector_values(Cache::DIFFUSION_COEFFS)[qp];

    const std::vector<libMesh::Real>& omega_dot = 
      cache.get_cached_vector_values(Cache::OMEGA_DOT)[qp];

    const libMesh::Number r = s_qpoint[qp](0);

    libMesh::Real jac = JxW[qp];

    if( this->_bc_handler->is_axisymmetric() )
      {
	jac *= r;
      }

    for(unsigned int s=0; s < this->_n_species; s++ )
      {
	libMesh::DenseSubVector<libMesh::Number> &Fs = 
	  *context.elem_subresiduals[this->_species_vars[s]]; // R_{s}

	const libMesh::Real term1 = -rho*(U*grad_w[s]) + omega_dot[s];
	const libMesh::Gradient term2 = -rho*D[s]*grad_w[s];

	for (unsigned int i=0; i != n_s_dofs; i++)
	  {
	    /*! \todo Need to add SCEBD term. */
	    Fs(i) += ( term1*s_phi[i][qp] + term2*s_grad_phi[i][qp] )*jac;

	    libmesh_assert( !libmesh_isnan(Fs(i)) );
	  }
      }

    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::assemble_momentum_time_deriv(libMesh::FEMContext& context, 
									  unsigned int qp,
									  const CachedValues& cache)
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.dof_indices_var[this->_u_var].size();

    // Check number of dofs is same for _u_var, v_var and w_var.
    libmesh_assert (n_u_dofs == context.dof_indices_var[this->_v_var].size());
    if (this->_dim == 3)
      libmesh_assert (n_u_dofs == context.dof_indices_var[this->_w_var].size());

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.element_fe_var[this->_u_var]->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.element_fe_var[this->_u_var]->get_phi();

    // The velocity shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.element_fe_var[this->_u_var]->get_dphi();

    const std::vector<libMesh::Point>& u_qpoint = 
      context.element_fe_var[this->_u_var]->get_xyz();

    libMesh::DenseSubVector<libMesh::Number> &Fu = *context.elem_subresiduals[this->_u_var]; // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = *context.elem_subresiduals[this->_v_var]; // R_{v}
    libMesh::DenseSubVector<libMesh::Number> &Fw = *context.elem_subresiduals[this->_w_var]; // R_{w}
    
    libMesh::Number rho = cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp];

    libMesh::Number u, v, w;
    u = cache.get_cached_values(Cache::X_VELOCITY)[qp];
    v = cache.get_cached_values(Cache::Y_VELOCITY)[qp];
    if (this->_dim == 3)
      w = cache.get_cached_values(Cache::Z_VELOCITY)[qp];

    libMesh::NumberVectorValue U(u,v);
    if (this->_dim == 3)
      U(2) = w;

    libMesh::Gradient grad_u = cache.get_cached_gradient_values(Cache::X_VELOCITY_GRAD)[qp];
    libMesh::Gradient grad_v = cache.get_cached_gradient_values(Cache::Y_VELOCITY_GRAD)[qp];

    libMesh::Gradient grad_w;
    if (this->_dim == 3)
      grad_w = cache.get_cached_gradient_values(Cache::Z_VELOCITY_GRAD)[qp];

    libMesh::Number divU = grad_u(0) + grad_v(1);
    if (this->_dim == 3)
      divU += grad_w(2);

    libMesh::NumberVectorValue grad_uT( grad_u(0), grad_v(0) ); 
    libMesh::NumberVectorValue grad_vT( grad_u(1), grad_v(1) );
    libMesh::NumberVectorValue grad_wT;
    if( this->_dim == 3 )
      {
	grad_uT(2) = grad_w(0);
	grad_vT(2) = grad_w(1);
	grad_wT = libMesh::NumberVectorValue( grad_u(2), grad_v(2), grad_w(2) );
      }

    libMesh::Number mu = cache.get_cached_values(Cache::MIXTURE_VISCOSITY)[qp];
    libMesh::Number p = cache.get_cached_values(Cache::PRESSURE)[qp];

    const libMesh::Number r = u_qpoint[qp](0);

    libMesh::Real jac = JxW[qp];

    if( this->_bc_handler->is_axisymmetric() )
      {
	divU += U(0)/r;
	jac *= r;
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
		   )*jac; 

	/*! \todo Would it be better to put this in its own DoF loop and do the if check once?*/
	if( this->_bc_handler->is_axisymmetric() )
	  {
	    Fu(i) += u_phi[i][qp]*( p/r - 2*mu*U(0)/(r*r) )*jac;
	  }

	libmesh_assert( !libmesh_isnan(Fu(i)) );

	Fv(i) += ( -rho*U*grad_v*u_phi[i][qp]                 // convection term
		   + p*u_gradphi[i][qp](1)                           // pressure term
		   - mu*(u_gradphi[i][qp]*grad_v + u_gradphi[i][qp]*grad_vT
			 - 2.0/3.0*divU*u_gradphi[i][qp](1) )    // diffusion term
		   + rho*this->_g(1)*u_phi[i][qp]                 // hydrostatic term
		   )*jac;

	
	libmesh_assert( !libmesh_isnan(Fv(i)) );

	if (this->_dim == 3)
	  {
	    Fw(i) += ( -rho*U*grad_w*u_phi[i][qp]                 // convection term
		       + p*u_gradphi[i][qp](2)                           // pressure term
		       - mu*(u_gradphi[i][qp]*grad_w + u_gradphi[i][qp]*grad_wT
			     - 2.0/3.0*divU*u_gradphi[i][qp](2) )    // diffusion term
		       + rho*this->_g(2)*u_phi[i][qp]                 // hydrostatic term
		       )*jac;

	    libmesh_assert( !libmesh_isnan(Fw(i)) );
	  }
      }
    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::assemble_energy_time_deriv( libMesh::FEMContext& context, 
									 unsigned int qp,
									 const CachedValues& cache)
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.dof_indices_var[this->_T_var].size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.element_fe_var[this->_T_var]->get_JxW();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.element_fe_var[this->_T_var]->get_phi();

    // The temperature shape functions gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.element_fe_var[this->_T_var]->get_dphi();

    // Physical location of the quadrature points
    const std::vector<libMesh::Point>& T_qpoint =
      context.element_fe_var[this->_T_var]->get_xyz();

    libMesh::DenseSubVector<libMesh::Number> &FT = *context.elem_subresiduals[this->_T_var]; // R_{T}

    libMesh::Number rho = cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp];

    libMesh::Number u, v, w;
    u = cache.get_cached_values(Cache::X_VELOCITY)[qp];
    v = cache.get_cached_values(Cache::Y_VELOCITY)[qp];
    if (this->_dim == 3)
      w = cache.get_cached_values(Cache::Z_VELOCITY)[qp];

    libMesh::NumberVectorValue U(u,v);
    if (this->_dim == 3)
      U(2) = w;

    libMesh::Number cp = 
      cache.get_cached_values(Cache::MIXTURE_SPECIFIC_HEAT_P)[qp];

    libMesh::Number k = 
      cache.get_cached_values(Cache::MIXTURE_THERMAL_CONDUCTIVITY)[qp];

    const libMesh::Gradient& grad_T = 
      cache.get_cached_gradient_values(Cache::TEMPERATURE_GRAD)[qp];

    const std::vector<libMesh::Real>& omega_dot = 
      cache.get_cached_vector_values(Cache::OMEGA_DOT)[qp];

    const std::vector<libMesh::Real>& h = 
      cache.get_cached_vector_values(Cache::SPECIES_ENTHALPY)[qp];

    libMesh::Real chem_term = 0.0;
    
    for(unsigned int s=0; s < this->_n_species; s++ )
      {
	chem_term += h[s]*omega_dot[s];
      }

    libmesh_assert( !libmesh_isnan(chem_term) );

    libMesh::Real jac = JxW[qp];

    if( this->_bc_handler->is_axisymmetric() )
      {
	const libMesh::Number r = T_qpoint[qp](0);
	jac *= r;
      }
    
    for (unsigned int i=0; i != n_T_dofs; i++)
      {
	FT(i) += ( ( -rho*cp*U*grad_T - chem_term )*T_phi[i][qp] // convection term + chemistry term
		     - k*grad_T*T_gradphi[i][qp]   /* diffusion term */   )*jac;

	libmesh_assert( !libmesh_isnan(FT(i)) );
      }

    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::compute_element_time_derivative_cache( const libMesh::FEMContext& context, 
										    CachedValues& cache ) const
  {
    const unsigned int n_qpoints = context.element_qrule->n_points();

    std::vector<libMesh::Real> u, v, w, T, p, p0;
    u.resize(n_qpoints);
    v.resize(n_qpoints);
    if( this->_dim > 2 )
      w.resize(n_qpoints);
    
    T.resize(n_qpoints);
    p.resize(n_qpoints);
    p0.resize(n_qpoints);

    std::vector<libMesh::Gradient> grad_u, grad_v, grad_w, grad_T;
    grad_u.resize(n_qpoints);
    grad_v.resize(n_qpoints);
    if( this->_dim > 2 )
      grad_w.resize(n_qpoints);
    
    grad_T.resize(n_qpoints);

    std::vector<std::vector<libMesh::Real> > mass_fractions;
    std::vector<std::vector<libMesh::Gradient> > grad_mass_fractions;
    mass_fractions.resize(n_qpoints);
    grad_mass_fractions.resize(n_qpoints);

    std::vector<libMesh::Real> M;
    M.resize(n_qpoints);

    std::vector<libMesh::Real> R;
    R.resize(n_qpoints);

    std::vector<libMesh::Real> rho;
    rho.resize(n_qpoints);

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
	u[qp] = context.interior_value(this->_u_var, qp);
	v[qp] = context.interior_value(this->_v_var, qp);

	grad_u[qp] = context.interior_gradient(this->_u_var, qp);
	grad_v[qp] = context.interior_gradient(this->_v_var, qp);
	if( this->_dim > 2 )
	  {
	    w[qp] = context.interior_value(this->_w_var, qp);
	    grad_w[qp] = context.interior_gradient(this->_w_var, qp);
	  }
	T[qp] = context.interior_value(this->_T_var, qp);
	grad_T[qp] = context.interior_gradient(this->_T_var, qp);

	p[qp] = context.interior_value(this->_p_var, qp);
	p0[qp] = this->get_p0_steady(context, qp);

	mass_fractions[qp].resize(this->_n_species);
	grad_mass_fractions[qp].resize(this->_n_species);

	for( unsigned int s = 0; s < this->_n_species; s++ )
	  {
	    /*! \todo Need to figure out something smarter for controling species
	              that go slightly negative. */
	    mass_fractions[qp][s] = std::max( context.interior_value(this->_species_vars[s],qp), 0.0 );
	    grad_mass_fractions[qp][s] = context.interior_gradient(this->_species_vars[s],qp);
	  }
	
	M[qp] = this->_gas_mixture.M( mass_fractions[qp] );

	R[qp] = this->_gas_mixture.R( mass_fractions[qp] );

	rho[qp] = this->rho( T[qp], p0[qp], mass_fractions[qp] );
      }
    
    cache.set_values(Cache::X_VELOCITY, u);
    cache.set_values(Cache::Y_VELOCITY, v);
    
    cache.set_gradient_values(Cache::X_VELOCITY_GRAD, grad_u);
    cache.set_gradient_values(Cache::Y_VELOCITY_GRAD, grad_v);
    
    if(this->_dim > 2)
      {
	cache.set_values(Cache::Z_VELOCITY, w);
	cache.set_gradient_values(Cache::Z_VELOCITY_GRAD, grad_w);
      }

    cache.set_values(Cache::TEMPERATURE, T);
    cache.set_gradient_values(Cache::TEMPERATURE_GRAD, grad_T);

    cache.set_values(Cache::PRESSURE, p);
    cache.set_values(Cache::THERMO_PRESSURE, p0);

    cache.set_vector_values(Cache::MASS_FRACTIONS, mass_fractions);
    cache.set_vector_gradient_values(Cache::MASS_FRACTIONS_GRAD, grad_mass_fractions);

    cache.set_values(Cache::MOLAR_MASS, M);

    cache.set_values(Cache::MIXTURE_GAS_CONSTANT, R);

    cache.set_values(Cache::MIXTURE_DENSITY, rho);

    /* These quantities must be computed after rho, T, mass_fractions, p0
       are set into the cache. */
    std::vector<libMesh::Real> mu;
    mu.resize(n_qpoints);

    std::vector<libMesh::Real> cp;
    cp.resize(n_qpoints);

    std::vector<libMesh::Real> k;
    k.resize(n_qpoints);

    std::vector<std::vector<libMesh::Real> > h;
    h.resize(n_qpoints);

    std::vector<std::vector<libMesh::Real> > h_RT_minus_s_R;
    h_RT_minus_s_R.resize(n_qpoints);

    std::vector<std::vector<libMesh::Real> > molar_densities;
    molar_densities.resize(n_qpoints);

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
	mu[qp] = this->_gas_mixture.mu(cache,qp);
	cp[qp] = this->_gas_mixture.cp(cache,qp);
	k[qp]  = this->_gas_mixture.k(cache,qp);

	h[qp].resize(this->_n_species);
	this->_gas_mixture.h( cache, qp, h[qp] );

	h_RT_minus_s_R[qp].resize(this->_n_species);
	this->_gas_mixture.h_RT_minus_s_R( cache, qp, h_RT_minus_s_R[qp] );

	molar_densities[qp].resize(this->_n_species);
	this->_gas_mixture.chem_mixture().molar_densities( rho[qp], mass_fractions[qp], 
							   molar_densities[qp] );
      }

    cache.set_values(Cache::MIXTURE_VISCOSITY, mu);
    cache.set_values(Cache::MIXTURE_SPECIFIC_HEAT_P, cp);
    cache.set_values(Cache::MIXTURE_THERMAL_CONDUCTIVITY, k);
    cache.set_vector_values(Cache::SPECIES_ENTHALPY, h);
    cache.set_vector_values(Cache::SPECIES_NORMALIZED_ENTHALPY_MINUS_NORMALIZED_ENTROPY, 
			    h_RT_minus_s_R);
    cache.set_vector_values(Cache::MOLAR_DENSITIES, molar_densities);

    /* Diffusion coefficients need rho, cp, k computed first.
       omega_dot may need h_RT_minus_s_R, molar_densities. */
    std::vector<std::vector<libMesh::Real> > D;
    D.resize(n_qpoints);

    std::vector<std::vector<libMesh::Real> > omega_dot;
    omega_dot.resize(n_qpoints);

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
	D[qp].resize(this->_n_species);
	this->_gas_mixture.D( cache, qp, D[qp] );

	omega_dot[qp].resize(this->_n_species);
	this->_gas_mixture.omega_dot( cache, qp, omega_dot[qp] );
      }

    cache.set_vector_values(Cache::DIFFUSION_COEFFS, D);
    cache.set_vector_values(Cache::OMEGA_DOT, omega_dot);

    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::compute_side_time_derivative_cache( const libMesh::FEMContext& context, 
										 CachedValues& cache ) const
  {
    const unsigned int n_qpoints = context.side_qrule->n_points();

    // Need for Catalytic Wall
    /*! \todo Add mechanism for checking if this side is a catalytic wall so we don't 
              compute these quantities unnecessarily. */
    std::vector<libMesh::Real> T, rho;
    T.resize(n_qpoints);
    rho.resize(n_qpoints);

    std::vector<std::vector<libMesh::Real> > mass_fractions;
    mass_fractions.resize(n_qpoints);

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
	T[qp] = context.side_value(this->_T_var, qp);

	mass_fractions[qp].resize(this->_n_species);
	for( unsigned int s = 0; s < this->_n_species; s++ )
	  {
	    /*! \todo Need to figure out something smarter for controling species
	              that go slightly negative. */
	    mass_fractions[qp][s] = std::max( context.side_value(this->_species_vars[s],qp), 0.0 );
	  }
	const libMesh::Real p0 = this->get_p0_steady_side(context, qp);

	rho[qp] = this->rho( T[qp], p0, mass_fractions[qp] );
      }

    cache.set_values(Cache::TEMPERATURE, T);
    cache.set_vector_values(Cache::MASS_FRACTIONS, mass_fractions);
    cache.set_values(Cache::MIXTURE_DENSITY, rho);

    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::compute_element_cache( const libMesh::FEMContext& context, 
								    const std::vector<libMesh::Point>& points,
								    CachedValues& cache ) const
  {
    if( cache.is_active(Cache::MIXTURE_DENSITY) )
      {
	std::vector<libMesh::Real> rho_values;
	rho_values.reserve( points.size() );

	std::vector<libMesh::Real> mass_fracs( this->_n_species );
	
	for( std::vector<libMesh::Point>::const_iterator point = points.begin();
	     point != points.end(); point++ )
	  {
	    libMesh::Real T = this->T(*point,context);
	    libMesh::Real p0 = this->get_p0_steady(context,*point);
	    this->mass_fractions( *point, context, mass_fracs );

	    rho_values.push_back(this->rho( T, p0, mass_fracs) );
	  }

	cache.set_values( Cache::MIXTURE_DENSITY, rho_values );
      }

    if( cache.is_active(Cache::SPECIES_VISCOSITY) )
      {
	libmesh_not_implemented();
      }

    if( cache.is_active(Cache::MOLE_FRACTIONS) )
      {
	std::vector<std::vector<libMesh::Real> > mole_fractions;
	mole_fractions.resize( points.size() );

	for( unsigned int i = 0; i < points.size(); i++ )
	  {
	    mole_fractions[i].resize( this->_n_species );
	  }

	std::vector<libMesh::Real> mass_fracs( this->_n_species );

	for( unsigned int p = 0; p < points.size(); p++ )
	  {
	    this->mass_fractions( points[p], context, mass_fracs );

	    libMesh::Real M = this->_gas_mixture.M(mass_fracs);

	    this->_gas_mixture.X(M,mass_fracs,mole_fractions[p]);
	  }

	cache.set_vector_values(Cache::MOLE_FRACTIONS, mole_fractions );
      }

    if( cache.is_active(Cache::SPECIES_ENTHALPY) )
      {
	{
	  std::vector<libMesh::Real> T;
	  T.resize( points.size() );

	  for( unsigned int p = 0; p < points.size(); p++ )
	    {
	      T[p] = this->T(points[p],context);
	    }

	  cache.set_values( Cache::TEMPERATURE, T );
	}

	std::vector<std::vector<libMesh::Real> > h;
	h.resize( points.size() );
	for( unsigned int p = 0; p < points.size(); p++ )
	  {
	    h[p].resize(this->_n_species);
	    this->_gas_mixture.h(cache, p, h[p]);
	  }

	cache.set_vector_values( Cache::SPECIES_ENTHALPY, h );
      }

    if( cache.is_active(Cache::OMEGA_DOT) )
      {
	{
	  std::vector<libMesh::Real> T, p0, rho, R;
	  T.resize( points.size() );
	  p0.resize( points.size() );
	  rho.resize( points.size() );
	  R.resize( points.size() );

	  std::vector<std::vector<libMesh::Real> > Y, molar_densities;
	  Y.resize( points.size() );
	  molar_densities.resize( points.size() );

	  for( unsigned int p = 0; p < points.size(); p++ )
	    {
	      T[p] = this->T(points[p],context);
	      p0[p] = this->get_p0_steady(context,points[p]);

	      Y[p].resize(this->_n_species);
	      this->mass_fractions( points[p], context, Y[p] );

	      rho[p] = this->rho( T[p], p0[p], Y[p] );

	      R[p] = this->_gas_mixture.R( Y[p] );

	      molar_densities[p].resize( this->_n_species );
	      this->_gas_mixture.chem_mixture().molar_densities( rho[p], 
								 Y[p], 
								 molar_densities[p] );
	    }

	  cache.set_values( Cache::TEMPERATURE, T );
	  cache.set_values( Cache::THERMO_PRESSURE, p0 );
	  cache.set_values( Cache::MIXTURE_DENSITY, rho );
	  cache.set_values( Cache::MIXTURE_GAS_CONSTANT, R );
	  cache.set_vector_values( Cache::MASS_FRACTIONS, Y );
	  cache.set_vector_values( Cache::MOLAR_DENSITIES, molar_densities );
	  
	  std::vector<std::vector<libMesh::Real> > h;
	  h.resize( points.size() );
	  for( unsigned int p = 0; p < points.size(); p++ )
	    {
	      h[p].resize(this->_n_species);
	      this->_gas_mixture.h_RT_minus_s_R( cache, p, h[p] );
	    }
	  cache.set_vector_values( Cache::SPECIES_NORMALIZED_ENTHALPY_MINUS_NORMALIZED_ENTROPY, h );

	}

	std::vector<std::vector<libMesh::Real> > omega_dot;
	omega_dot.resize( points.size() );

	for( unsigned int p = 0; p < points.size(); p++ )
	  {
	    omega_dot[p].resize(this->_n_species);
	    this->_gas_mixture.omega_dot( cache, p, omega_dot[p] );
	  }

	cache.set_vector_values(Cache::OMEGA_DOT, omega_dot );
      }

    return;
  }

  // Instantiate
  template class ReactingLowMachNavierStokes< IdealGasMixture<CEAThermodynamics,ConstantTransport,Kinetics> >;
#ifdef GRINS_HAVE_CANTERA
  template class ReactingLowMachNavierStokes< IdealGasMixture<CanteraThermodynamics,CanteraTransport,CanteraKinetics> >;
  template class ReactingLowMachNavierStokes< IdealGasMixture<CanteraThermodynamics,ConstantTransport,CanteraKinetics> >;
  template class ReactingLowMachNavierStokes< IdealGasMixture<CEAThermodynamics,ConstantTransport,CanteraKinetics> >;
#endif

} // namespace GRINS
