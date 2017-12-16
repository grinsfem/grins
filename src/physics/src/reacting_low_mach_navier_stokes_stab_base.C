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
#include "grins/reacting_low_mach_navier_stokes_stab_base.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/physical_constants.h"

namespace GRINS
{

  template<typename Mixture, typename Evaluator>
  ReactingLowMachNavierStokesStabilizationBase<Mixture,Evaluator>::ReactingLowMachNavierStokesStabilizationBase
  ( const std::string& physics_name,const GetPot& input,std::unique_ptr<Mixture> & gas_mix )
    : ReactingLowMachNavierStokesBase<Mixture>(physics_name,input,gas_mix),
    _stab_helper( physics_name+"StabHelper", input )
  {}

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokesStabilizationBase<Mixture,Evaluator>::init_context( AssemblyContext& context )
  {
    // First call base class
    ReactingLowMachNavierStokesAbstract::init_context(context);

    // We need pressure derivatives
    context.get_element_fe(this->_press_var.p())->get_dphi();

    // We also need second derivatives, so initialize those.
    context.get_element_fe(this->_flow_vars.u())->get_d2phi();
    context.get_element_fe(this->_temp_vars.T())->get_d2phi();
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokesStabilizationBase<Mixture,Evaluator>::compute_res_steady( AssemblyContext& context,
                                                                                            unsigned int qp,
                                                                                            libMesh::Real& RP_s,
                                                                                            libMesh::RealGradient& RM_s,
                                                                                            libMesh::Real& RE_s,
                                                                                            std::vector<libMesh::Real>& Rs_s )
  {
    Rs_s.resize(this->n_species(),0.0);

    // Grab r-coordinate for axisymmetric terms
    // We're assuming all variables are using the same quadrature rule
    libMesh::Real r = (context.get_element_fe(this->_flow_vars.u())->get_xyz())[qp](0);

    libMesh::RealGradient grad_p = context.interior_gradient(this->_press_var.p(), qp);

    libMesh::RealGradient grad_u = context.interior_gradient(this->_flow_vars.u(), qp);

    libMesh::RealGradient U( context.interior_value(this->_flow_vars.u(), qp) );
    libMesh::Real divU = grad_u(0);

    if( this->_is_axisymmetric )
      divU += U(0)/r;

    if(this->_flow_vars.dim() > 1)
      {
        U(1) = context.interior_value(this->_flow_vars.v(), qp);
        divU += (context.interior_gradient(this->_flow_vars.v(), qp))(1);
      }
    if(this->_flow_vars.dim() == 3)
      {
        U(2) = context.interior_value(this->_flow_vars.w(), qp);
        divU += (context.interior_gradient(this->_flow_vars.w(), qp))(2);
      }

    // We don't add axisymmetric terms here since we don't directly use hess_{u,v}
    // axisymmetric terms are built into divGradU, etc. functions below
    libMesh::RealTensor hess_u = context.interior_hessian(this->_flow_vars.u(), qp);
    libMesh::RealTensor hess_v = context.interior_hessian(this->_flow_vars.v(), qp);

    libMesh::Real T = context.interior_value(this->_temp_vars.T(), qp);

    libMesh::Gradient grad_T = context.interior_gradient(this->_temp_vars.T(), qp);
    libMesh::Tensor hess_T = context.interior_hessian(this->_temp_vars.T(), qp);

    libMesh::Real hess_T_term = hess_T(0,0) + hess_T(1,1);
#if LIBMESH_DIM > 2
    hess_T_term += hess_T(2,2);
#endif
    // Add axisymmetric terms, if needed
    if( this->_is_axisymmetric )
      hess_T_term += grad_T(0)/r;

    std::vector<libMesh::Real> ws(this->n_species());
    std::vector<libMesh::RealGradient> grad_ws(this->n_species());
    std::vector<libMesh::RealTensor> hess_ws(this->n_species());
    for(unsigned int s=0; s < this->_n_species; s++ )
      {
        ws[s] = context.interior_value(this->_species_vars.species(s), qp);
        grad_ws[s] = context.interior_gradient(this->_species_vars.species(s), qp);
        hess_ws[s] = context.interior_hessian(this->_species_vars.species(s), qp);
      }

    Evaluator gas_evaluator( *(this->_gas_mixture) );
    const libMesh::Real R_mix = gas_evaluator.R_mix(ws);
    const libMesh::Real p0 = this->get_p0_steady(context,qp);
    libMesh::Real rho = this->rho(T, p0, R_mix );
    libMesh::Real cp = gas_evaluator.cp(T,p0,ws);
    libMesh::Real M = gas_evaluator.M_mix( ws );

    std::vector<libMesh::Real> D( this->n_species() );
    libMesh::Real mu, k;

    gas_evaluator.mu_and_k_and_D( T, rho, cp, ws, mu, k, D );


    // grad_rho = drho_dT*gradT + \sum_s drho_dws*grad_ws
    const libMesh::Real drho_dT = -p0/(R_mix*T*T);
    libMesh::RealGradient grad_rho = drho_dT*grad_T;
    for(unsigned int s=0; s < this->_n_species; s++ )
      {
        libMesh::Real Ms = gas_evaluator.M(s);
        libMesh::Real R_uni = Constants::R_universal/1000.0; /* J/kmol-K --> J/mol-K */

        // drho_dws = -p0/(T*R_mix*R_mix)*dR_dws
        // dR_dws = R_uni*d_dws(1/M)
        // d_dws(1/M) = d_dws(\sum_s w_s/Ms) =  1/Ms
        const libMesh::Real drho_dws = -p0/(R_mix*R_mix*T)*R_uni/Ms;
        grad_rho += drho_dws*grad_ws[s];
      }

    libMesh::RealGradient rhoUdotGradU;
    libMesh::RealGradient divGradU;
    libMesh::RealGradient divGradUT;
    libMesh::RealGradient divdivU;

    if( this->_flow_vars.dim() == 1 )
      {
        rhoUdotGradU = rho*_stab_helper.UdotGradU( U, grad_u );

        divGradU  = _stab_helper.div_GradU( hess_u );
        divGradUT = _stab_helper.div_GradU_T( hess_u );
        divdivU   = _stab_helper.div_divU_I( hess_u );
      }
    else if( this->_flow_vars.dim() == 2 )
      {
        libMesh::RealGradient grad_v = context.interior_gradient(this->_flow_vars.v(), qp);

        rhoUdotGradU = rho*_stab_helper.UdotGradU( U, grad_u, grad_v );

        // Call axisymmetric versions if we are doing an axisymmetric run
        if( this->_is_axisymmetric )
          {
            divGradU  = _stab_helper.div_GradU_axi( r, U, grad_u, grad_v, hess_u, hess_v );
            divGradUT = _stab_helper.div_GradU_T_axi( r, U, grad_u, hess_u, hess_v );
            divdivU   = _stab_helper.div_divU_I_axi( r, U, grad_u, hess_u, hess_v );
          }
        else
          {
            divGradU  = _stab_helper.div_GradU( hess_u, hess_v );
            divGradUT = _stab_helper.div_GradU_T( hess_u, hess_v );
            divdivU   = _stab_helper.div_divU_I( hess_u, hess_v );
          }
      }
    else
      {
        libMesh::RealGradient grad_v = context.interior_gradient(this->_flow_vars.v(), qp);
        libMesh::RealTensor hess_v = context.interior_hessian(this->_flow_vars.v(), qp);

        libMesh::RealGradient grad_w = context.interior_gradient(this->_flow_vars.w(), qp);
        libMesh::RealTensor hess_w = context.interior_hessian(this->_flow_vars.w(), qp);

        rhoUdotGradU = rho*_stab_helper.UdotGradU( U, grad_u, grad_v, grad_w );

        divGradU  = _stab_helper.div_GradU( hess_u, hess_v, hess_w );
        divGradUT = _stab_helper.div_GradU_T( hess_u, hess_v, hess_w );
        divdivU   = _stab_helper.div_divU_I( hess_u, hess_v, hess_w );
      }



    // Terms if we have vicosity derivatives w.r.t. temp.
    /*
      if( this->_mu.deriv(T) != 0.0 )
      {
      libMesh::Gradient gradTgradu( grad_T*grad_u, grad_T*grad_v );

      libMesh::Gradient gradTgraduT( grad_T(0)*grad_u(0) + grad_T(1)*grad_u(1),
      grad_T(0)*grad_v(0) + grad_T(1)*grad_v(1) );

      libMesh::Real divU = grad_u(0) + grad_v(1);

      libMesh::Gradient gradTdivU( grad_T(0)*divU, grad_T(1)*divU );

      if(this->_flow_vars.dim() == 3)
      {
      libMesh::Gradient grad_w = context.interior_gradient(this->_flow_vars.w(), qp);

      gradTgradu(2) = grad_T*grad_w;

      gradTgraduT(0) += grad_T(2)*grad_u(2);
      gradTgraduT(1) += grad_T(2)*grad_v(2);
      gradTgraduT(2) = grad_T(0)*grad_w(0) + grad_T(1)*grad_w(1) + grad_T(2)*grad_w(2);

      divU += grad_w(2);
      gradTdivU(0) += grad_T(0)*grad_w(2);
      gradTdivU(1) += grad_T(1)*grad_w(2);
      gradTdivU(2) += grad_T(2)*divU;
      }

      divT += this->_mu.deriv(T)*( gradTgradu + gradTgraduT - 2.0/3.0*gradTdivU );
      }
    */

    // Axisymmetric terms already built in
    libMesh::RealGradient div_stress = mu*(divGradU + divGradUT - 2.0/3.0*divdivU);

    std::vector<libMesh::Real> omega_dot(this->n_species());
    gas_evaluator.omega_dot(T,rho,ws,omega_dot);

    libMesh::Real chem_term = 0.0;
    libMesh::Gradient mass_term(0.0,0.0,0.0);
    for(unsigned int s=0; s < this->_n_species; s++ )
      {
        // Start accumulating chemistry term for energy residual
        libMesh::Real h_s=gas_evaluator.h_s(T,s);
        chem_term += h_s*omega_dot[s];

        /* Accumulate mass term for continuity residual
           mass_term = grad_M/M */
        mass_term += grad_ws[s]/this->_gas_mixture->M(s);

        libMesh::Real hess_s_term = hess_ws[s](0,0) + hess_ws[s](1,1);
#if LIBMESH_DIM > 2
        hess_s_term += hess_ws[s](2,2);
#endif
        // Add axisymmetric terms, if needed
        if( this->_is_axisymmetric )
          hess_s_term += grad_ws[s](0)/r;

        // Species residual
        /*! \todo Still missing derivative of species diffusion coefficient.
          rho*grad_D[s]*grad_ws[s] */
        Rs_s[s] = rho*U*grad_ws[s] - rho*D[s]*hess_s_term - grad_rho*D[s]*grad_ws[s]
          - omega_dot[s];
      }
    mass_term *= M;

    // Continuity residual
    RP_s = divU - (U*grad_T)/T - U*mass_term;

    // Momentum residual
    RM_s = rhoUdotGradU + grad_p - div_stress - rho*(this->_g);

    // Energy residual
    // - this->_k.deriv(T)*(grad_T*grad_T)
    RE_s = rho*U*cp*grad_T  - k*(hess_T_term) + chem_term;

    return;
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokesStabilizationBase<Mixture,Evaluator>::compute_res_transient( AssemblyContext& context,
                                                                                               unsigned int qp,
                                                                                               libMesh::Real& RP_t,
                                                                                               libMesh::RealGradient& RM_t,
                                                                                               libMesh::Real& RE_t,
                                                                                               std::vector<libMesh::Real>& Rs_t )
  {
    libMesh::Real T = context.interior_value( this->_temp_vars.T(), qp );

    std::vector<libMesh::Real> ws(this->n_species());
    for(unsigned int s=0; s < this->_n_species; s++ )
      {
        ws[s] = context.interior_value(this->_species_vars.species(s), qp);
      }

    Evaluator gas_evaluator( *(this->_gas_mixture) );
    const libMesh::Real R_mix = gas_evaluator.R_mix(ws);
    const libMesh::Real p0 = this->get_p0_transient(context,qp);
    const libMesh::Real rho = this->rho(T, p0, R_mix);
    const libMesh::Real cp = gas_evaluator.cp(T,p0,ws);
    const libMesh::Real M = gas_evaluator.M_mix( ws );

    // M_dot = -M^2 \sum_s w_dot[s]/Ms
    libMesh::Real M_dot = 0.0;
    std::vector<libMesh::Real> ws_dot(this->n_species());
    for(unsigned int s=0; s < this->n_species(); s++)
      {
        context.interior_rate(this->_species_vars.species(s), qp, ws_dot[s]);

        // Start accumulating M_dot
        M_dot += ws_dot[s]/this->_gas_mixture->M(s);
      }
    libMesh::Real M_dot_over_M = M_dot*(-M);

    libMesh::RealGradient u_dot;
    context.interior_rate(this->_flow_vars.u(), qp, u_dot(0));
    if(this->_flow_vars.dim() > 1)
      context.interior_rate(this->_flow_vars.v(), qp, u_dot(1));
    if(this->_flow_vars.dim() == 3)
      context.interior_rate(this->_flow_vars.w(), qp, u_dot(2));

    libMesh::Real T_dot;
    context.interior_rate(this->_temp_vars.T(), qp, T_dot);

    RP_t = -T_dot/T + M_dot_over_M;
    RM_t = rho*u_dot;
    RE_t = rho*cp*T_dot;
    for(unsigned int s=0; s < this->n_species(); s++)
      {
        Rs_t[s] = rho*ws_dot[s];
      }

    return;
  }

} // namespace GRINS
