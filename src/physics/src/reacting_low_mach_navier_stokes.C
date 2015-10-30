//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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
#include "grins/postprocessed_quantities.h"

// libMesh
#include "libmesh/quadrature.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  template<typename Mixture, typename Evaluator>
  ReactingLowMachNavierStokes<Mixture,Evaluator>::ReactingLowMachNavierStokes(const PhysicsName& physics_name, const GetPot& input)
    : ReactingLowMachNavierStokesBase<Mixture,Evaluator>(physics_name,input),
    _p_pinning(input,physics_name),
    _rho_index(0),
    _mu_index(0),
    _k_index(0),
    _cp_index(0)
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
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::auxiliary_init( MultiphysicsSystem& system )
  {
    if( _pin_pressure )
      _p_pinning.check_pin_location(system.get_mesh());
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::register_postprocessing_vars( const GetPot& input,
                                                                                     PostProcessedQuantities<libMesh::Real>& postprocessing )
  {
    std::string section = "Physics/"+reacting_low_mach_navier_stokes+"/output_vars";

    if( input.have_variable(section) )
      {
        unsigned int n_vars = input.vector_variable_size(section);

        for( unsigned int v = 0; v < n_vars; v++ )
          {
            std::string name = input(section,"DIE!",v);

            if( name == std::string("rho") )
              {
                this->_rho_index = postprocessing.register_quantity( name );
              }
            else if( name == std::string("mu") )
              {
                this->_mu_index = postprocessing.register_quantity( name );
              }
            else if( name == std::string("k") )
              {
                this->_k_index = postprocessing.register_quantity( name );
              }
            else if( name == std::string("cp") )
              {
                this->_cp_index = postprocessing.register_quantity( name );
              }
            else if( name == std::string("mole_fractions") )
              {
                this->_mole_fractions_index.resize(this->n_species());

                for( unsigned int s = 0; s < this->n_species(); s++ )
                  {
                    this->_mole_fractions_index[s] = postprocessing.register_quantity( "X_"+this->_gas_mixture.species_name(s) );
                  }
              }
            else if( name == std::string("h_s") )
              {
                this->_h_s_index.resize(this->n_species());

                for( unsigned int s = 0; s < this->n_species(); s++ )
                  {
                    this->_h_s_index[s] = postprocessing.register_quantity( "h_"+this->_gas_mixture.species_name(s) );
                  }
              }
            else if( name == std::string("omega_dot") )
              {
                this->_omega_dot_index.resize(this->n_species());

                for( unsigned int s = 0; s < this->n_species(); s++ )
                  {
                    this->_omega_dot_index[s] = postprocessing.register_quantity( "omega_dot_"+this->_gas_mixture.species_name(s) );
                  }

                std::cout << "omega_dot size = " << _omega_dot_index.size() << std::endl;
              }
            else
              {
                std::cerr << "Error: Invalue output_vars value for "+reacting_low_mach_navier_stokes << std::endl
                          << "       Found " << name << std::endl
                          << "       Acceptable values are: rho" << std::endl
                          << "                              mu" << std::endl
                          << "                              k" << std::endl
                          << "                              cp" << std::endl
                          << "                              mole_fractions" << std::endl
                          << "                              omega_dot" << std::endl;
                libmesh_error();
              }
          }
      }

    return;
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::init_context( AssemblyContext& context )
  {
    // First call base class
    ReactingLowMachNavierStokesBase<Mixture,Evaluator>::init_context(context);

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
    if( compute_jacobian )
      {
        libmesh_not_implemented();
      }
    // Convenience
    const VariableIndex s0_var = this->_species_vars[0];

    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.get_dof_indices(this->_p_var).size();
    const unsigned int n_s_dofs = context.get_dof_indices(s0_var).size();
    const unsigned int n_u_dofs = context.get_dof_indices(this->_u_var).size();
    const unsigned int n_T_dofs = context.get_dof_indices(this->_T_var).size();

    // Check number of dofs is same for _u_var, v_var and w_var.
    libmesh_assert (n_u_dofs == context.get_dof_indices(this->_v_var).size());
    if (this->_dim == 3)
      libmesh_assert (n_u_dofs == context.get_dof_indices(this->_w_var).size());

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real>& JxW =
      context.get_element_fe(this->_u_var)->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& p_phi =
      context.get_element_fe(this->_p_var)->get_phi();

    // The species shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& s_phi = context.get_element_fe(s0_var)->get_phi();

    // The species shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::Gradient> >& s_grad_phi = context.get_element_fe(s0_var)->get_dphi();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_u_var)->get_phi();

    // The velocity shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_u_var)->get_dphi();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_T_var)->get_phi();

    // The temperature shape functions gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(this->_T_var)->get_dphi();

    const std::vector<libMesh::Point>& u_qpoint =
      context.get_element_fe(this->_u_var)->get_xyz();

    libMesh::DenseSubVector<libMesh::Number>& Fp = context.get_elem_residual(this->_p_var); // R_{p}

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_u_var); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_v_var); // R_{v}
    libMesh::DenseSubVector<libMesh::Number> &Fw = context.get_elem_residual(this->_w_var); // R_{w}

    libMesh::DenseSubVector<libMesh::Number> &FT = context.get_elem_residual(this->_T_var); // R_{T}

    unsigned int n_qpoints = context.get_element_qrule().n_points();
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Number rho = cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp];

        libMesh::Number u, v, T;
        u = cache.get_cached_values(Cache::X_VELOCITY)[qp];
        v = cache.get_cached_values(Cache::Y_VELOCITY)[qp];

        T = cache.get_cached_values(Cache::TEMPERATURE)[qp];

        const libMesh::Gradient& grad_T =
          cache.get_cached_gradient_values(Cache::TEMPERATURE_GRAD)[qp];

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

        const std::vector<libMesh::Real>& D =
          cache.get_cached_vector_values(Cache::DIFFUSION_COEFFS)[qp];

        libMesh::Number cp =
          cache.get_cached_values(Cache::MIXTURE_SPECIFIC_HEAT_P)[qp];

        libMesh::Number k =
          cache.get_cached_values(Cache::MIXTURE_THERMAL_CONDUCTIVITY)[qp];

        const std::vector<libMesh::Real>& omega_dot =
          cache.get_cached_vector_values(Cache::OMEGA_DOT)[qp];

        const std::vector<libMesh::Real>& h =
          cache.get_cached_vector_values(Cache::SPECIES_ENTHALPY)[qp];

        const libMesh::Number r = u_qpoint[qp](0);

        libMesh::Real jac = JxW[qp];

        if( this->_is_axisymmetric )
          {
            divU += U(0)/r;
            jac *= r;
          }

        libMesh::Real M = cache.get_cached_values(Cache::MOLAR_MASS)[qp];

        std::vector<libMesh::Gradient> grad_ws = cache.get_cached_vector_gradient_values(Cache::MASS_FRACTIONS_GRAD)[qp];
        libmesh_assert_equal_to( grad_ws.size(), this->_n_species );

        // Continuity Residual
        libMesh::Gradient mass_term(0.0,0.0,0.0);
        for(unsigned int s=0; s < this->_n_species; s++ )
          {
            mass_term += grad_ws[s]/this->_gas_mixture.M(s);
          }
        mass_term *= M;

        for (unsigned int i=0; i != n_p_dofs; i++)
          {
            Fp(i) += (-U*(mass_term + grad_T/T) + divU)*jac*p_phi[i][qp];
          }

        // Species residuals
        for(unsigned int s=0; s < this->_n_species; s++ )
          {
            libMesh::DenseSubVector<libMesh::Number> &Fs =
              context.get_elem_residual(this->_species_vars[s]); // R_{s}

            const libMesh::Real term1 = -rho*(U*grad_ws[s]) + omega_dot[s];
            const libMesh::Gradient term2 = -rho*D[s]*grad_ws[s];

            for (unsigned int i=0; i != n_s_dofs; i++)
              {
                /*! \todo Need to add SCEBD term. */
                Fs(i) += ( term1*s_phi[i][qp] + term2*s_grad_phi[i][qp] )*jac;
              }
          }

        // Momentum residuals
        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) += ( -rho*U*grad_u*u_phi[i][qp]
                       + p*u_gradphi[i][qp](0)
                       - mu*(u_gradphi[i][qp]*grad_u + u_gradphi[i][qp]*grad_uT
                             - 2.0/3.0*divU*u_gradphi[i][qp](0) )
                       + rho*this->_g(0)*u_phi[i][qp]
                       )*jac;

            /*! \todo Would it be better to put this in its own DoF loop
                      and do the if check once?*/
            if( this->_is_axisymmetric )
              {
                Fu(i) += u_phi[i][qp]*( p/r - 2*mu*U(0)/(r*r) - 2.0/3.0*mu*divU/r )*jac;
              }

            Fv(i) += ( -rho*U*grad_v*u_phi[i][qp]
                       + p*u_gradphi[i][qp](1)
                       - mu*(u_gradphi[i][qp]*grad_v + u_gradphi[i][qp]*grad_vT
                             - 2.0/3.0*divU*u_gradphi[i][qp](1) )
                       + rho*this->_g(1)*u_phi[i][qp]
                       )*jac;

            if (this->_dim == 3)
              {
                Fw(i) += ( -rho*U*grad_w*u_phi[i][qp]
                           + p*u_gradphi[i][qp](2)
                           - mu*(u_gradphi[i][qp]*grad_w + u_gradphi[i][qp]*grad_wT
                                 - 2.0/3.0*divU*u_gradphi[i][qp](2) )
                           + rho*this->_g(2)*u_phi[i][qp]
                           )*jac;
              }
          }

        // Energy residual
        libMesh::Real chem_term = 0.0;
        for(unsigned int s=0; s < this->_n_species; s++ )
          {
            chem_term += h[s]*omega_dot[s];
          }

        for (unsigned int i=0; i != n_T_dofs; i++)
          {
            FT(i) += ( ( -rho*cp*U*grad_T - chem_term )*T_phi[i][qp]
                       - k*grad_T*T_gradphi[i][qp]  )*jac;
          }

      } // quadrature loop

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
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::element_constraint( bool compute_jacobian,
                                                                           AssemblyContext& context,
                                                                           CachedValues& /* cache */ )
  {
    // Pin p = p_value at p_point
    if( this->_pin_pressure )
      {
	_p_pinning.pin_value( context, compute_jacobian, this->_p_var );
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
        cp[qp] = gas_evaluator.cp(cache,qp);

        D_s[qp].resize(this->_n_species);

        gas_evaluator.mu_and_k_and_D( T[qp], rho[qp], cp[qp], mass_fractions[qp],
                                      mu[qp], k[qp], D_s[qp] );

	h_s[qp].resize(this->_n_species);
	gas_evaluator.h_s( cache, qp, h_s[qp] );

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
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::compute_postprocessed_quantity( unsigned int quantity_index,
                                                                                       const AssemblyContext& context,
                                                                                       const libMesh::Point& point,
                                                                                       libMesh::Real& value )
  {
    Evaluator gas_evaluator( this->_gas_mixture );

    if( quantity_index == this->_rho_index )
      {
        std::vector<libMesh::Real> Y( this->_n_species );
        libMesh::Real T = this->T(point,context);
        libMesh::Real p0 = this->get_p0_steady(context,point);
        this->mass_fractions( point, context, Y );

        value = this->rho( T, p0, gas_evaluator.R_mix(Y) );
      }
    else if( quantity_index == this->_mu_index )
      {
        std::vector<libMesh::Real> Y( this->_n_species );
        libMesh::Real T = this->T(point,context);
        this->mass_fractions( point, context, Y );

        value = gas_evaluator.mu( T, Y );
      }
    else if( quantity_index == this->_k_index )
      {
        std::vector<libMesh::Real> Y( this->_n_species );
        libMesh::Real T = this->T(point,context);
        this->mass_fractions( point, context, Y );

        value = gas_evaluator.k( T, Y );
      }
    else if( quantity_index == this->_cp_index )
      {
        std::vector<libMesh::Real> Y( this->_n_species );
        libMesh::Real T = this->T(point,context);
        this->mass_fractions( point, context, Y );

        value = gas_evaluator.cp( T, Y );
      }
    // Now check for species dependent stuff
    else
      {
        if( !this->_h_s_index.empty() )
          {
            libmesh_assert_equal_to( _h_s_index.size(), this->n_species() );

            for( unsigned int s = 0; s < this->n_species(); s++ )
              {
                if( quantity_index == this->_h_s_index[s] )
                  {
                    libMesh::Real T = this->T(point,context);

                    value = gas_evaluator.h_s( T, s );
                    return;
                  }
              }
          }

        if( !this->_mole_fractions_index.empty() )
          {
            libmesh_assert_equal_to( _mole_fractions_index.size(), this->n_species() );

            for( unsigned int s = 0; s < this->n_species(); s++ )
              {
                if( quantity_index == this->_mole_fractions_index[s] )
                  {
                    std::vector<libMesh::Real> Y( this->_n_species );
                    this->mass_fractions( point, context, Y );

                    libMesh::Real M = gas_evaluator.M_mix(Y);

                    value = gas_evaluator.X( s, M, Y[s] );
                    return;
                  }
              }
          }

        if( !this->_omega_dot_index.empty() )
          {
            libmesh_assert_equal_to( _omega_dot_index.size(), this->n_species() );

            for( unsigned int s = 0; s < this->n_species(); s++ )
              {
                if( quantity_index == this->_omega_dot_index[s] )
                  {
                    std::vector<libMesh::Real> Y( this->n_species() );
                    this->mass_fractions( point, context, Y );

                    libMesh::Real T = this->T(point,context);

                    libMesh::Real p0 = this->get_p0_steady(context,point);

                    libMesh::Real rho = this->rho( T, p0, gas_evaluator.R_mix(Y) );

                    std::vector<libMesh::Real> omega_dot( this->n_species() );
                    gas_evaluator.omega_dot( T, rho, Y, omega_dot );

                    value = omega_dot[s];
                    return;
                  }
              }
          }

      } // if/else quantity_index

    return;
  }

} // namespace GRINS
