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
#include "grins/reacting_low_mach_navier_stokes.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/cached_quantities_enum.h"
#include "grins/generic_ic_handler.h"
#include "grins/reacting_low_mach_navier_stokes_bc_handling.h"

// libMesh
#include "libmesh/quadrature.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  template<typename Mixture, typename Evaluator>
  ReactingLowMachNavierStokes<Mixture,Evaluator>::ReactingLowMachNavierStokes(const PhysicsName& physics_name, const GetPot& input)
    : ReactingLowMachNavierStokesBase(physics_name,input),
      _gas_mixture(input),
      _p_pinning(input,physics_name)
  {
    this->read_input_options(input);

    // This is deleted in the base class
    this->_bc_handler = new ReactingLowMachNavierStokesBCHandling<typename Mixture::ChemistryParent>( physics_name, input,
                                                                                                      this->_gas_mixture.chemistry() );

    if( this->_bc_handler->is_axisymmetric() )
      {
        this->_is_axisymmetric = true;
      }

    this->_ic_handler = new GenericICHandler( physics_name, input );

    return;
  }

  template<typename Mixture, typename Evaluator>
  ReactingLowMachNavierStokes<Mixture,Evaluator>::~ReactingLowMachNavierStokes()
  {
    return;
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::read_input_options( const GetPot& input )
  {
    // Other quantities read in base class

    // Read pressure pinning information
    this->_pin_pressure = input("Physics/"+reacting_low_mach_navier_stokes+"/pin_pressure", false );
  
    return;
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::init_context( AssemblyContext& context )
  {
    // First call base class
    GRINS::ReactingLowMachNavierStokesBase::init_context(context);

    // We also need the side shape functions, etc.
    context.get_side_fe(this->_u_var)->get_JxW();
    context.get_side_fe(this->_u_var)->get_phi();
    context.get_side_fe(this->_u_var)->get_dphi();
    context.get_side_fe(this->_u_var)->get_xyz();

    context.get_side_fe(this->_T_var)->get_JxW();
    context.get_side_fe(this->_T_var)->get_phi();
    context.get_side_fe(this->_T_var)->get_dphi();
    context.get_side_fe(this->_T_var)->get_xyz();

    return;
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::element_time_derivative( bool compute_jacobian,
                                                                                AssemblyContext& context,
                                                                                CachedValues& cache )
  {
    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	this->assemble_mass_time_deriv(context, qp, cache);
	this->assemble_species_time_deriv(context, qp, cache);
	this->assemble_momentum_time_deriv(context, qp, cache);
	this->assemble_energy_time_deriv(context, qp, cache);
      }

    return;
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::side_time_derivative( bool compute_jacobian,
                                                                             AssemblyContext& context,
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

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::side_constraint( bool compute_jacobian,
                                                                        AssemblyContext& context,
                                                                        CachedValues& /* cache */ )
  {
    // Pin p = p_value at p_point
    if( _pin_pressure )
      {
	_p_pinning.pin_value( context, compute_jacobian, _p_var );
      }

    return;
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::mass_residual( bool /*compute_jacobian*/,
                                                                      AssemblyContext& /*context*/,
                                                                      CachedValues& /*cache*/ )
  {
    libmesh_not_implemented();

    return;
  }


  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::assemble_mass_time_deriv( AssemblyContext& context, 
                                                                                 unsigned int qp,
                                                                                 const CachedValues& cache )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.get_dof_indices(this->_p_var).size();
    
    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real>& JxW =
      context.get_element_fe(this->_u_var)->get_JxW();
    
    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& p_phi =
      context.get_element_fe(this->_p_var)->get_phi();
    
    const std::vector<libMesh::Point>& u_qpoint = 
      context.get_element_fe(this->_u_var)->get_xyz();

    libMesh::DenseSubVector<libMesh::Number>& Fp = context.get_elem_residual(this->_p_var); // R_{p}

    libMesh::Number u, v, T;
    u = cache.get_cached_values(Cache::X_VELOCITY)[qp];
    v = cache.get_cached_values(Cache::Y_VELOCITY)[qp];
    
    T = cache.get_cached_values(Cache::TEMPERATURE)[qp];

    libMesh::NumberVectorValue U(u,v);
    if (this->_dim == 3)
      U(2) = cache.get_cached_values(Cache::Z_VELOCITY)[qp]; // w

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

    if( this->_is_axisymmetric )
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

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::assemble_species_time_deriv(AssemblyContext& context, 
                                                                                   unsigned int qp,
                                                                                   const CachedValues& cache)
  {
    // Convenience
    const VariableIndex s0_var = this->_species_vars[0];
    
    /* The number of local degrees of freedom in each species variable.
       We assume the same number of dofs for each species */
    const unsigned int n_s_dofs = context.get_dof_indices(s0_var).size();
    
    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW = context.get_element_fe(s0_var)->get_JxW();
    
    // The species shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& s_phi = context.get_element_fe(s0_var)->get_phi();

    // The species shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::Gradient> >& s_grad_phi = context.get_element_fe(s0_var)->get_dphi();

    const std::vector<libMesh::Point>& s_qpoint = 
      context.get_element_fe(this->_species_vars[0])->get_xyz();

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

    if( this->_is_axisymmetric )
      {
	jac *= r;
      }

    for(unsigned int s=0; s < this->_n_species; s++ )
      {
	libMesh::DenseSubVector<libMesh::Number> &Fs = 
	  context.get_elem_residual(this->_species_vars[s]); // R_{s}

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

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::assemble_momentum_time_deriv(AssemblyContext& context, 
									  unsigned int qp,
									  const CachedValues& cache)
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.get_dof_indices(this->_u_var).size();

    // Check number of dofs is same for _u_var, v_var and w_var.
    libmesh_assert (n_u_dofs == context.get_dof_indices(this->_v_var).size());
    if (this->_dim == 3)
      libmesh_assert (n_u_dofs == context.get_dof_indices(this->_w_var).size());

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_u_var)->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_u_var)->get_phi();

    // The velocity shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_u_var)->get_dphi();

    const std::vector<libMesh::Point>& u_qpoint = 
      context.get_element_fe(this->_u_var)->get_xyz();

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_u_var); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_v_var); // R_{v}
    libMesh::DenseSubVector<libMesh::Number> &Fw = context.get_elem_residual(this->_w_var); // R_{w}
    
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

    if( this->_is_axisymmetric )
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
	if( this->_is_axisymmetric )
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

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::assemble_energy_time_deriv( AssemblyContext& context, 
                                                                                   unsigned int qp,
                                                                                   const CachedValues& cache)
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.get_dof_indices(this->_T_var).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_T_var)->get_JxW();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_T_var)->get_phi();

    // The temperature shape functions gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(this->_T_var)->get_dphi();

    // Physical location of the quadrature points
    const std::vector<libMesh::Point>& T_qpoint =
      context.get_element_fe(this->_T_var)->get_xyz();

    libMesh::DenseSubVector<libMesh::Number> &FT = context.get_elem_residual(this->_T_var); // R_{T}

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

    if( this->_is_axisymmetric )
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

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::compute_element_time_derivative_cache( const AssemblyContext& context, 
                                                                                              CachedValues& cache )
  {
    Evaluator gas_evaluator( this->_gas_mixture );

    const unsigned int n_qpoints = context.get_element_qrule().n_points();

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
	
	M[qp] = gas_evaluator.M_mix( mass_fractions[qp] );

	R[qp] = gas_evaluator.R_mix( mass_fractions[qp] );

	rho[qp] = this->rho( T[qp], p0[qp], R[qp] );
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

    /* These quantities must be computed after T, mass_fractions, p0
       are set into the cache. */
    std::vector<libMesh::Real> mu;
    mu.resize(n_qpoints);

    std::vector<libMesh::Real> cp;
    cp.resize(n_qpoints);

    std::vector<libMesh::Real> k;
    k.resize(n_qpoints);

    std::vector<std::vector<libMesh::Real> > h_s;
    h_s.resize(n_qpoints);

    std::vector<std::vector<libMesh::Real> > D_s;
    D_s.resize(n_qpoints);

    std::vector<std::vector<libMesh::Real> > omega_dot_s;
    omega_dot_s.resize(n_qpoints);

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
	gas_evaluator.mu_and_k(cache,qp,mu[qp],k[qp]);
	cp[qp] = gas_evaluator.cp(cache,qp);

	h_s[qp].resize(this->_n_species);
	gas_evaluator.h_s( cache, qp, h_s[qp] );

	D_s[qp].resize(this->_n_species);
	gas_evaluator.D( cache, qp, D_s[qp] );

	omega_dot_s[qp].resize(this->_n_species);
	gas_evaluator.omega_dot( cache, qp, omega_dot_s[qp] );
      }

    cache.set_values(Cache::MIXTURE_VISCOSITY, mu);
    cache.set_values(Cache::MIXTURE_SPECIFIC_HEAT_P, cp);
    cache.set_values(Cache::MIXTURE_THERMAL_CONDUCTIVITY, k);
    cache.set_vector_values(Cache::SPECIES_ENTHALPY, h_s);
    cache.set_vector_values(Cache::DIFFUSION_COEFFS, D_s);
    cache.set_vector_values(Cache::OMEGA_DOT, omega_dot_s);

    return;
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::compute_side_time_derivative_cache( const AssemblyContext& context, 
                                                                                           CachedValues& cache )
  {
    Evaluator gas_evaluator( this->_gas_mixture );

    const unsigned int n_qpoints = context.get_side_qrule().n_points();

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

	rho[qp] = this->rho( T[qp], p0, gas_evaluator.R_mix(mass_fractions[qp]) );
      }

    cache.set_values(Cache::TEMPERATURE, T);
    cache.set_vector_values(Cache::MASS_FRACTIONS, mass_fractions);
    cache.set_values(Cache::MIXTURE_DENSITY, rho);

    return;
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::compute_element_cache( const AssemblyContext& context, 
                                                                              const std::vector<libMesh::Point>& points,
                                                                              CachedValues& cache )
  {
    Evaluator gas_evaluator( this->_gas_mixture );

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

	    rho_values.push_back(this->rho( T, p0, gas_evaluator.R_mix(mass_fracs) ) );
	  }

	cache.set_values( Cache::MIXTURE_DENSITY, rho_values );
      }

    if( cache.is_active(Cache::SPECIES_VISCOSITY) )
      {
	libmesh_not_implemented();
      }

    if( cache.is_active(Cache::MIXTURE_VISCOSITY) )
      {
	std::vector<libMesh::Real> mu_values;
	mu_values.reserve( points.size() );

        std::vector<libMesh::Real> T;
        T.resize( points.size() );

	std::vector<std::vector<libMesh::Real> >  mass_fracs;
	mass_fracs.resize( points.size() );

        for( unsigned int i = 0; i < points.size(); i++ )
	  {
            mass_fracs[i].resize( this->_n_species );
	  }

	for( unsigned int p = 0; p < points.size(); p++ )
	  {
	    T[p] = this->T(points[p],context);

            this->mass_fractions( points[p], context, mass_fracs[p] );
          }
        
        cache.set_values( Cache::TEMPERATURE, T );
        cache.set_vector_values(Cache::MASS_FRACTIONS, mass_fracs );

        for( unsigned int p = 0; p < points.size(); p++ )
	  {
            mu_values.push_back( gas_evaluator.mu( cache, p ) );
	  }

	cache.set_values( Cache::MIXTURE_VISCOSITY, mu_values );
      }

    if( cache.is_active(Cache::SPECIES_THERMAL_CONDUCTIVITY) )
      {
	libmesh_not_implemented();
      }

    if( cache.is_active(Cache::MIXTURE_THERMAL_CONDUCTIVITY) )
      {
	std::vector<libMesh::Real> k_values;
	k_values.reserve( points.size() );

        std::vector<libMesh::Real> T;
        T.resize( points.size() );

	std::vector<std::vector<libMesh::Real> >  mass_fracs;
	mass_fracs.resize( points.size() );

        for( unsigned int i = 0; i < points.size(); i++ )
	  {
            mass_fracs[i].resize( this->_n_species );
	  }

	for( unsigned int p = 0; p < points.size(); p++ )
	  {
	    T[p] = this->T(points[p],context);

            this->mass_fractions( points[p], context, mass_fracs[p] );
          }
        
        cache.set_values( Cache::TEMPERATURE, T );
        cache.set_vector_values(Cache::MASS_FRACTIONS, mass_fracs );

        for( unsigned int p = 0; p < points.size(); p++ )
	  {
	    k_values.push_back( gas_evaluator.k( cache, p ) );
	  }

	cache.set_values( Cache::MIXTURE_THERMAL_CONDUCTIVITY, k_values );
      }

    if( cache.is_active(Cache::SPECIES_SPECIFIC_HEAT_P) )
      {
	libmesh_not_implemented();
      }

    if( cache.is_active(Cache::MIXTURE_SPECIFIC_HEAT_P) )
      {
	std::vector<libMesh::Real> cp_values;
	cp_values.reserve( points.size() );

        std::vector<libMesh::Real> T;
        T.resize( points.size() );

	std::vector<std::vector<libMesh::Real> >  mass_fracs;
	mass_fracs.resize( points.size() );

        for( unsigned int i = 0; i < points.size(); i++ )
	  {
            mass_fracs[i].resize( this->_n_species );
	  }

	for( unsigned int p = 0; p < points.size(); p++ )
	  {
	    T[p] = this->T(points[p],context);

            this->mass_fractions( points[p], context, mass_fracs[p] );
          }
        
        cache.set_values( Cache::TEMPERATURE, T );
        cache.set_vector_values(Cache::MASS_FRACTIONS, mass_fracs );

        for( unsigned int p = 0; p < points.size(); p++ )
	  {
	    cp_values.push_back( gas_evaluator.cp( cache, p ) );
	  }

	cache.set_values( Cache::MIXTURE_SPECIFIC_HEAT_P, cp_values );
      }

    if( cache.is_active(Cache::SPECIES_SPECIFIC_HEAT_V) )
      {
	libmesh_not_implemented();
      }

    if( cache.is_active(Cache::MIXTURE_SPECIFIC_HEAT_V) )
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

	    libMesh::Real M = gas_evaluator.M_mix(mass_fracs);

	    gas_evaluator.X(M,mass_fracs,mole_fractions[p]);
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
	    gas_evaluator.h_s(cache, p, h[p]);
	  }

	cache.set_vector_values( Cache::SPECIES_ENTHALPY, h );
      }

    if( cache.is_active(Cache::OMEGA_DOT) )
      {
	{
	  std::vector<libMesh::Real> T, R, rho, p0;
	  T.resize( points.size() );
	  R.resize( points.size() );
          rho.resize( points.size() );
          p0.resize( points.size() );

	  std::vector<std::vector<libMesh::Real> > Y;
	  Y.resize( points.size() );

	  for( unsigned int p = 0; p < points.size(); p++ )
	    {
	      T[p] = this->T(points[p],context);

	      Y[p].resize(this->_n_species);
	      this->mass_fractions( points[p], context, Y[p] );

              R[p] = gas_evaluator.R_mix( Y[p] );

              p0[p] = this->get_p0_steady(context, points[p] );

              rho[p] = this->rho( T[p], p0[p], R[p] );
	    }

	  cache.set_values( Cache::TEMPERATURE, T );
	  cache.set_values( Cache::MIXTURE_GAS_CONSTANT, R);
          cache.set_values(Cache::MIXTURE_DENSITY, rho);
          cache.set_values(Cache::THERMO_PRESSURE, p0);
	  cache.set_vector_values( Cache::MASS_FRACTIONS, Y );
	}

	std::vector<std::vector<libMesh::Real> > omega_dot;
	omega_dot.resize( points.size() );

	for( unsigned int p = 0; p < points.size(); p++ )
	  {
	    omega_dot[p].resize(this->_n_species);
	    gas_evaluator.omega_dot( cache, p, omega_dot[p] );
	  }

	cache.set_vector_values(Cache::OMEGA_DOT, omega_dot );
      }

    return;
  }
  
  template<typename Mixture, typename Evaluator>
  libMesh::Real ReactingLowMachNavierStokes<Mixture,Evaluator>::cp_mix( const libMesh::Real T,
                                                                        const std::vector<libMesh::Real>& Y )
  {
    Evaluator gas_evaluator( this->_gas_mixture );
    
    return gas_evaluator.cp( T, Y );
  }

  template<typename Mixture, typename Evaluator>
  libMesh::Real ReactingLowMachNavierStokes<Mixture,Evaluator>::mu( const libMesh::Real T,
                                                                    const std::vector<libMesh::Real>& Y )
  {
    Evaluator gas_evaluator( this->_gas_mixture );
    
    return gas_evaluator.mu( T, Y );
  }

  template<typename Mixture, typename Evaluator>
  libMesh::Real ReactingLowMachNavierStokes<Mixture,Evaluator>::k( const libMesh::Real T,
                                                                   const std::vector<libMesh::Real>& Y )
  {
    Evaluator gas_evaluator( this->_gas_mixture );
    
    return gas_evaluator.k( T, Y );
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::D( const libMesh::Real rho,
                                                          const libMesh::Real cp,
                                                          const libMesh::Real k,
                                                          std::vector<libMesh::Real>& D )
  {
    Evaluator gas_evaluator( this->_gas_mixture );
    
    return gas_evaluator.D( rho, cp, k, D );
  }

} // namespace GRINS
