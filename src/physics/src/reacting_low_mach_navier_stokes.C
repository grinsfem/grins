//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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
#include "grins/postprocessed_quantities.h"

// libMesh
#include "libmesh/quadrature.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  template<typename Mixture, typename Evaluator>
  ReactingLowMachNavierStokes<Mixture,Evaluator>::ReactingLowMachNavierStokes(const PhysicsName& physics_name,
                                                                              const GetPot& input,
                                                                              std::unique_ptr<Mixture> & gas_mix)
    : ReactingLowMachNavierStokesBase<Mixture>(physics_name,input,gas_mix),
    _p_pinning(input,physics_name),
    _rho_index(0),
    _mu_index(0),
    _k_index(0),
    _cp_index(0)
  {
    this->_pin_pressure = input("Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/pin_pressure", false );

    this->_ic_handler = new GenericICHandler( physics_name, input );
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
    std::string section = "Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/output_vars";

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
                    this->_mole_fractions_index[s] = postprocessing.register_quantity( "X_"+this->_gas_mixture->species_name(s) );
                  }
              }
            else if( name == std::string("h_s") )
              {
                this->_h_s_index.resize(this->n_species());

                for( unsigned int s = 0; s < this->n_species(); s++ )
                  {
                    this->_h_s_index[s] = postprocessing.register_quantity( "h_"+this->_gas_mixture->species_name(s) );
                  }
              }
            else if( name == std::string("omega_dot") )
              {
                this->_omega_dot_index.resize(this->n_species());

                for( unsigned int s = 0; s < this->n_species(); s++ )
                  this->_omega_dot_index[s] = postprocessing.register_quantity( "omega_dot_"+this->_gas_mixture->species_name(s) );
              }
            else if( name == std::string("D_s") )
              {
                this->_Ds_index.resize(this->n_species());

                for( unsigned int s = 0; s < this->n_species(); s++ )
                  this->_Ds_index[s] = postprocessing.register_quantity( "D_"+this->_gas_mixture->species_name(s) );
              }
            else
              {
                std::cerr << "Error: Invalue output_vars value for "+PhysicsNaming::reacting_low_mach_navier_stokes() << std::endl
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
    ReactingLowMachNavierStokesAbstract::init_context(context);

    // We also need the side shape functions, etc.
    context.get_side_fe(this->_flow_vars.u())->get_JxW();
    context.get_side_fe(this->_flow_vars.u())->get_phi();
    context.get_side_fe(this->_flow_vars.u())->get_dphi();
    context.get_side_fe(this->_flow_vars.u())->get_xyz();

    context.get_side_fe(this->_temp_vars.T())->get_JxW();
    context.get_side_fe(this->_temp_vars.T())->get_phi();
    context.get_side_fe(this->_temp_vars.T())->get_dphi();
    context.get_side_fe(this->_temp_vars.T())->get_xyz();
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::element_time_derivative
  ( bool compute_jacobian, AssemblyContext & context )
  {
    if( compute_jacobian )
      libmesh_not_implemented();

    const CachedValues & cache = context.get_cached_values();

    // Convenience
    const VariableIndex s0_var = this->_species_vars.species(0);

    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.get_dof_indices(this->_press_var.p()).size();
    const unsigned int n_s_dofs = context.get_dof_indices(s0_var).size();
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();

    // Check number of dofs is same for _flow_vars.u(), v_var and w_var.
    if (this->_flow_vars.dim() > 1)
      libmesh_assert (n_u_dofs == context.get_dof_indices(this->_flow_vars.v()).size());

    if (this->_flow_vars.dim() == 3)
      libmesh_assert (n_u_dofs == context.get_dof_indices(this->_flow_vars.w()).size());

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real>& JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& p_phi =
      context.get_element_fe(this->_press_var.p())->get_phi();

    // The species shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& s_phi = context.get_element_fe(s0_var)->get_phi();

    // The species shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::Gradient> >& s_grad_phi = context.get_element_fe(s0_var)->get_dphi();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    // The velocity shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_flow_vars.u())->get_dphi();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_temp_vars.T())->get_phi();

    // The temperature shape functions gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(this->_temp_vars.T())->get_dphi();

    const std::vector<libMesh::Point>& u_qpoint =
      context.get_element_fe(this->_flow_vars.u())->get_xyz();

    libMesh::DenseSubVector<libMesh::Number>& Fp = context.get_elem_residual(this->_press_var.p()); // R_{p}

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u()); // R_{u}

    libMesh::DenseSubVector<libMesh::Number>* Fv = NULL;
    if( this->_flow_vars.dim() > 1 )
      Fv = &context.get_elem_residual(this->_flow_vars.v()); // R_{v}

    libMesh::DenseSubVector<libMesh::Number>* Fw = NULL;
    if( this->_flow_vars.dim() == 3 )
      Fw = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}

    libMesh::DenseSubVector<libMesh::Number> &FT = context.get_elem_residual(this->_temp_vars.T()); // R_{T}

    unsigned int n_qpoints = context.get_element_qrule().n_points();
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Number rho = cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp];

        libMesh::Number u, T;
        u = cache.get_cached_values(Cache::X_VELOCITY)[qp];

        T = cache.get_cached_values(Cache::TEMPERATURE)[qp];

        const libMesh::Gradient& grad_T =
          cache.get_cached_gradient_values(Cache::TEMPERATURE_GRAD)[qp];

        libMesh::NumberVectorValue U(u);
        if (this->_flow_vars.dim() > 1)
          U(1) = cache.get_cached_values(Cache::Y_VELOCITY)[qp];
        if (this->_flow_vars.dim() == 3)
          U(2) = cache.get_cached_values(Cache::Z_VELOCITY)[qp]; // w

        libMesh::Gradient grad_u = cache.get_cached_gradient_values(Cache::X_VELOCITY_GRAD)[qp];
        libMesh::Gradient grad_v;

        if (this->_flow_vars.dim() > 1)
          grad_v = cache.get_cached_gradient_values(Cache::Y_VELOCITY_GRAD)[qp];

        libMesh::Gradient grad_w;
        if (this->_flow_vars.dim() == 3)
          grad_w = cache.get_cached_gradient_values(Cache::Z_VELOCITY_GRAD)[qp];

        libMesh::Number divU = grad_u(0);
        if (this->_flow_vars.dim() > 1)
          divU += grad_v(1);
        if (this->_flow_vars.dim() == 3)
          divU += grad_w(2);

        libMesh::NumberVectorValue grad_uT( grad_u(0) );
        libMesh::NumberVectorValue grad_vT;
        if( this->_flow_vars.dim() > 1 )
          {
            grad_uT(1) = grad_v(0);
            grad_vT = libMesh::NumberVectorValue( grad_u(1), grad_v(1) );
          }
        libMesh::NumberVectorValue grad_wT;
        if( this->_flow_vars.dim() == 3 )
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

        if(Physics::is_axisymmetric())
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
            mass_term += grad_ws[s]/this->_gas_mixture->M(s);
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
              context.get_elem_residual(this->_species_vars.species(s)); // R_{s}

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
            if(Physics::is_axisymmetric())
              {
                Fu(i) += u_phi[i][qp]*( p/r - 2*mu*U(0)/(r*r) - 2.0/3.0*mu*divU/r )*jac;
              }

            if (this->_flow_vars.dim() > 1)
              {
                (*Fv)(i) += ( -rho*U*grad_v*u_phi[i][qp]
                              + p*u_gradphi[i][qp](1)
                              - mu*(u_gradphi[i][qp]*grad_v + u_gradphi[i][qp]*grad_vT
                                    - 2.0/3.0*divU*u_gradphi[i][qp](1) )
                              + rho*this->_g(1)*u_phi[i][qp]
                              )*jac;
              }

            if (this->_flow_vars.dim() == 3)
              {
                (*Fw)(i) += ( -rho*U*grad_w*u_phi[i][qp]
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
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::element_constraint
  ( bool compute_jacobian,
    AssemblyContext & context )
  {
    // Pin p = p_value at p_point
    if( this->_pin_pressure )
      {
        _p_pinning.pin_value( context, compute_jacobian, this->_press_var.p() );
      }

    return;
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::mass_residual
  ( bool compute_jacobian, AssemblyContext & context )
  {
    const unsigned int n_p_dofs = context.get_dof_indices(this->_press_var.p()).size();
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();
    const VariableIndex s0_var = this->_species_vars.species(0);
    const unsigned int n_s_dofs = context.get_dof_indices(s0_var).size();

    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    const std::vector<std::vector<libMesh::Real> >& p_phi =
      context.get_element_fe(this->_press_var.p())->get_phi();

    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_temp_vars.T())->get_phi();

    const std::vector<std::vector<libMesh::Real> >& s_phi =
      context.get_element_fe(s0_var)->get_phi();


    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<libMesh::Real> &F_p = context.get_elem_residual(this->_press_var.p());

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<libMesh::Real> &F_u = context.get_elem_residual(this->_flow_vars.u());

    libMesh::DenseSubVector<libMesh::Real>* F_v = NULL;
    if( this->_flow_vars.dim() > 1 )
      F_v  = &context.get_elem_residual(this->_flow_vars.v());

    libMesh::DenseSubVector<libMesh::Real>* F_w = NULL;
    if( this->_flow_vars.dim() == 3 )
      F_w  = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}

    libMesh::DenseSubVector<libMesh::Real> &F_T = context.get_elem_residual(this->_temp_vars.T());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    const std::vector<libMesh::Point>& u_qpoint =
      context.get_element_fe(this->_flow_vars.u())->get_xyz();

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
        libMesh::Real u_dot, v_dot = 0.0, w_dot = 0.0;
        context.interior_rate(this->_flow_vars.u(), qp, u_dot);

        if( this->_flow_vars.dim() > 1 )
          context.interior_rate(this->_flow_vars.v(), qp, v_dot);

        if( this->_flow_vars.dim() == 3 )
          context.interior_rate(this->_flow_vars.w(), qp, w_dot);

        libMesh::Real T_dot;
        context.interior_rate(this->_temp_vars.T(), qp, T_dot);

        libMesh::Real T = context.interior_value(this->_temp_vars.T(), qp);

        std::vector<libMesh::Real> ws(this->n_species());
        for(unsigned int s=0; s < this->_n_species; s++ )
          ws[s] = context.interior_value(this->_species_vars.species(s), qp);

        Evaluator gas_evaluator( *(this->_gas_mixture) );
        const libMesh::Real R_mix = gas_evaluator.R_mix(ws);
        const libMesh::Real p0 = this->get_p0_steady(context,qp);
        const libMesh::Real rho = this->rho(T, p0, R_mix);
        const libMesh::Real cp = gas_evaluator.cp(T,p0,ws);
        const libMesh::Real M = gas_evaluator.M_mix( ws );

        libMesh::Real jac = JxW[qp];
        const libMesh::Number r = u_qpoint[qp](0);

        if( this->_is_axisymmetric )
          jac *= r;

        // M_dot = -M^2 \sum_s w_dot[s]/Ms
        libMesh::Real M_dot = 0.0;

        // Species residual
        for(unsigned int s=0; s < this->n_species(); s++)
          {
            libMesh::DenseSubVector<libMesh::Number> &F_s =
              context.get_elem_residual(this->_species_vars.species(s));

            libMesh::Real ws_dot;
            context.interior_rate(this->_species_vars.species(s), qp, ws_dot);

            for (unsigned int i = 0; i != n_s_dofs; ++i)
              F_s(i) -= rho*ws_dot*s_phi[i][qp]*jac;

            // Start accumulating M_dot
            M_dot += ws_dot/this->_gas_mixture->M(s);
          }

        // Continuity residual
        // M_dot = -M^2 \sum_s w_dot[s]/Ms
        libMesh::Real M_dot_over_M = M_dot*(-M);

        for (unsigned int i = 0; i != n_p_dofs; ++i)
          F_p(i) -= (T_dot/T-M_dot_over_M)*p_phi[i][qp]*jac;

        // Momentum residual
        for (unsigned int i = 0; i != n_u_dofs; ++i)
          {
            F_u(i) -= rho*u_dot*u_phi[i][qp]*jac;

            if( this->_flow_vars.dim() > 1 )
              (*F_v)(i) -= rho*v_dot*u_phi[i][qp]*jac;

            if( this->_flow_vars.dim() == 3 )
              (*F_w)(i) -= rho*w_dot*u_phi[i][qp]*jac;
          }

        // Energy residual
        for (unsigned int i = 0; i != n_T_dofs; ++i)
          F_T(i) -= rho*cp*T_dot*T_phi[i][qp]*jac;

        if( compute_jacobian )
          libmesh_not_implemented();

      }
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::compute_element_time_derivative_cache
  ( AssemblyContext & context )
  {
    CachedValues & cache = context.get_cached_values();

    Evaluator gas_evaluator( *(this->_gas_mixture) );

    const unsigned int n_qpoints = context.get_element_qrule().n_points();

    std::vector<libMesh::Real> u, v, w, T, p, p0;
    u.resize(n_qpoints);
    v.resize(n_qpoints);
    if( this->_flow_vars.dim() > 2 )
      w.resize(n_qpoints);

    T.resize(n_qpoints);
    p.resize(n_qpoints);
    p0.resize(n_qpoints);

    std::vector<libMesh::Gradient> grad_u, grad_v, grad_w, grad_T;
    grad_u.resize(n_qpoints);
    grad_v.resize(n_qpoints);
    if( this->_flow_vars.dim() > 2 )
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

    std::vector<libMesh::Real> cp;
    cp.resize(n_qpoints);

    std::vector<libMesh::Real> mu;
    mu.resize(n_qpoints);

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
        u[qp] = context.interior_value(this->_flow_vars.u(), qp);
        v[qp] = context.interior_value(this->_flow_vars.v(), qp);

        grad_u[qp] = context.interior_gradient(this->_flow_vars.u(), qp);
        grad_v[qp] = context.interior_gradient(this->_flow_vars.v(), qp);
        if( this->_flow_vars.dim() > 2 )
          {
            w[qp] = context.interior_value(this->_flow_vars.w(), qp);
            grad_w[qp] = context.interior_gradient(this->_flow_vars.w(), qp);
          }

        T[qp] = context.interior_value(this->_temp_vars.T(), qp);
        grad_T[qp] = context.interior_gradient(this->_temp_vars.T(), qp);

        p[qp] = context.interior_value(this->_press_var.p(), qp);
        p0[qp] = this->get_p0_steady(context, qp);

        mass_fractions[qp].resize(this->_n_species);
        grad_mass_fractions[qp].resize(this->_n_species);
        h_s[qp].resize(this->_n_species);

        for( unsigned int s = 0; s < this->_n_species; s++ )
          {
            /*! \todo Need to figure out something smarter for controling species
              that go slightly negative. */
            mass_fractions[qp][s] = std::max( context.interior_value(this->_species_vars.species(s),qp), 0.0 );
            grad_mass_fractions[qp][s] = context.interior_gradient(this->_species_vars.species(s),qp);
            h_s[qp][s] = gas_evaluator.h_s( T[qp], s );
          }

        M[qp] = gas_evaluator.M_mix( mass_fractions[qp] );

        R[qp] = gas_evaluator.R_mix( mass_fractions[qp] );

        rho[qp] = this->rho( T[qp], p0[qp], R[qp] );

        cp[qp] = gas_evaluator.cp(T[qp], p0[qp], mass_fractions[qp]);

        D_s[qp].resize(this->_n_species);

        gas_evaluator.mu_and_k_and_D( T[qp], rho[qp], cp[qp], mass_fractions[qp],
                                      mu[qp], k[qp], D_s[qp] );

        omega_dot_s[qp].resize(this->_n_species);
        gas_evaluator.omega_dot( T[qp], rho[qp], mass_fractions[qp], omega_dot_s[qp] );
      }

    cache.set_values(Cache::X_VELOCITY, u);
    cache.set_values(Cache::Y_VELOCITY, v);
    cache.set_gradient_values(Cache::X_VELOCITY_GRAD, grad_u);
    cache.set_gradient_values(Cache::Y_VELOCITY_GRAD, grad_v);

    if(this->_flow_vars.dim() > 2)
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
    cache.set_values(Cache::MIXTURE_SPECIFIC_HEAT_P, cp);
    cache.set_values(Cache::MIXTURE_VISCOSITY, mu);
    cache.set_values(Cache::MIXTURE_THERMAL_CONDUCTIVITY, k);
    cache.set_vector_values(Cache::DIFFUSION_COEFFS, D_s);
    cache.set_vector_values(Cache::SPECIES_ENTHALPY, h_s);
    cache.set_vector_values(Cache::OMEGA_DOT, omega_dot_s);
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokes<Mixture,Evaluator>::compute_postprocessed_quantity( unsigned int quantity_index,
                                                                                       const AssemblyContext& context,
                                                                                       const libMesh::Point& point,
                                                                                       libMesh::Real& value )
  {
    Evaluator gas_evaluator( *(this->_gas_mixture) );

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
        libMesh::Real p0 = this->get_p0_steady(context,point);

        value = gas_evaluator.mu( T, p0, Y );
      }
    else if( quantity_index == this->_k_index )
      {
        std::vector<libMesh::Real> Y( this->_n_species );

        libMesh::Real T = this->T(point,context);
        this->mass_fractions( point, context, Y );
        libMesh::Real p0 = this->get_p0_steady(context,point);

        libMesh::Real cp = gas_evaluator.cp( T, p0, Y );

        libMesh::Real rho = this->rho( T, p0, gas_evaluator.R_mix(Y) );

        libMesh::Real mu, k;
        std::vector<libMesh::Real> D( this->_n_species );

        gas_evaluator.mu_and_k_and_D( T, rho, cp, Y, mu, k, D );

        value = k;
        return;
      }
    else if( quantity_index == this->_cp_index )
      {
        std::vector<libMesh::Real> Y( this->_n_species );
        libMesh::Real T = this->T(point,context);
        this->mass_fractions( point, context, Y );
        libMesh::Real p0 = this->get_p0_steady(context,point);

        value = gas_evaluator.cp( T, p0, Y );
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

        if( !this->_Ds_index.empty() )
          {
            libmesh_assert_equal_to( _Ds_index.size(), this->n_species() );

            for( unsigned int s = 0; s < this->n_species(); s++ )
              {
                if( quantity_index == this->_Ds_index[s] )
                  {
                    std::vector<libMesh::Real> Y( this->_n_species );

                    libMesh::Real T = this->T(point,context);
                    this->mass_fractions( point, context, Y );
                    libMesh::Real p0 = this->get_p0_steady(context,point);

                    libMesh::Real cp = gas_evaluator.cp( T, p0, Y );

                    libMesh::Real rho = this->rho( T, p0, gas_evaluator.R_mix(Y) );

                    libMesh::Real mu, k;
                    std::vector<libMesh::Real> D( this->_n_species );

                    gas_evaluator.mu_and_k_and_D( T, rho, cp, Y, mu, k, D );

                    value = D[s];
                    return;
                  }
              }
          }

      } // if/else quantity_index

    return;
  }

} // namespace GRINS
