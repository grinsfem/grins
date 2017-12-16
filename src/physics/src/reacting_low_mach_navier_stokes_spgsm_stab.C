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
#include "grins/reacting_low_mach_navier_stokes_spgsm_stab.h"

// GRINS
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{
  template<typename Mixture, typename Evaluator>
  ReactingLowMachNavierStokesSPGSMStabilization<Mixture,Evaluator>::ReactingLowMachNavierStokesSPGSMStabilization
  ( const GRINS::PhysicsName& physics_name, const GetPot& input,std::unique_ptr<Mixture> & gas_mix )
    : ReactingLowMachNavierStokesStabilizationBase<Mixture,Evaluator>(physics_name,input,gas_mix)
  {}

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokesSPGSMStabilization<Mixture,Evaluator>::element_time_derivative
  ( bool compute_jacobian,
    AssemblyContext & context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.get_dof_indices(this->_press_var.p()).size();
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();
    const VariableIndex s0_var = this->_species_vars.species(0);
    const unsigned int n_s_dofs = context.get_dof_indices(s0_var).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
      context.get_element_fe(this->_press_var.p())->get_dphi();

    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_flow_vars.u())->get_dphi();

    // The temperature shape functions gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(this->_temp_vars.T())->get_dphi();

    const std::vector<std::vector<libMesh::Gradient> >& s_gradphi = context.get_element_fe(s0_var)->get_dphi();

    // We're assuming the quadrature rule is the same for all variables
    const std::vector<libMesh::Point>& u_qpoint =
      context.get_element_fe(this->_flow_vars.u())->get_xyz();

    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(this->_press_var.p()); // R_{p}

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u()); // R_{u}

    libMesh::DenseSubVector<libMesh::Number>* Fv = NULL;
    if( this->_flow_vars.dim() > 1 )
      Fv  = &context.get_elem_residual(this->_flow_vars.v()); // R_{v}

    libMesh::DenseSubVector<libMesh::Number>* Fw = NULL;
    if( this->_flow_vars.dim() == 3 )
      Fw  = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}

    libMesh::DenseSubVector<libMesh::Number> &FT = context.get_elem_residual(this->_temp_vars.T()); // R_{T}

    unsigned int n_qpoints = context.get_element_qrule().n_points();


    libMesh::FEBase* u_fe = context.get_element_fe(this->_flow_vars.u());
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Real T = context.interior_value( this->_temp_vars.T(), qp );

        libMesh::RealGradient U( context.interior_value( this->_flow_vars.u(), qp ) );
        if( this->_flow_vars.dim() > 1 )
          U(1) = context.interior_value( this->_flow_vars.v(), qp );
        if( this->_flow_vars.dim() == 3 )
          U(2) = context.interior_value( this->_flow_vars.w(), qp );

        std::vector<libMesh::Real> ws(this->n_species());
        for(unsigned int s=0; s < this->_n_species; s++ )
          {
            ws[s] = context.fixed_interior_value(this->_species_vars.species(s), qp);
          }

        Evaluator gas_evaluator( *(this->_gas_mixture) );
        const libMesh::Real R_mix = gas_evaluator.R_mix(ws);
        const libMesh::Real p0 = this->get_p0_steady(context,qp);
        libMesh::Real rho = this->rho(T, p0, R_mix);

        const libMesh::Real cp = gas_evaluator.cp(T,p0,ws);

        std::vector<libMesh::Real> D( this->n_species() );
        libMesh::Real mu, k;

        gas_evaluator.mu_and_k_and_D( T, rho, cp, ws, mu, k, D );

        libMesh::RealGradient g = this->_stab_helper.compute_g( u_fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( u_fe, context, qp );

        // Taus
        libMesh::Real tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, rho, U, mu, this->_is_steady );

        libMesh::Real tau_C = this->_stab_helper.compute_tau_continuity( tau_M, g );

        libMesh::Real tau_E = this->_stab_helper.compute_tau_energy( context, qp, g, G, rho, U, k, cp, this->_is_steady );



        // Strong form residuals
        libMesh::RealGradient RM_s = 0.0;
        libMesh::Real RC_s = 0.0;
        libMesh::Real RE_s = 0.0;
        std::vector<libMesh::Real> Rs_s;

        this->compute_res_steady( context, qp, RC_s, RM_s, RE_s, Rs_s );

        const libMesh::Number r = u_qpoint[qp](0);

        libMesh::Real jac = JxW[qp];

        if( this->_is_axisymmetric )
          jac *= r;

        // Pressure PSPG term
        for (unsigned int i=0; i != n_p_dofs; i++)
          Fp(i) += tau_M*RM_s*p_dphi[i][qp]*jac;

        // Momentum SUPG + div-div terms
        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) += ( - tau_C*RC_s*u_gradphi[i][qp](0)
                       - tau_M*RM_s(0)*rho*U*u_gradphi[i][qp]  )*jac;

            if( this->_is_axisymmetric )
              Fu(i) += (-tau_C*RC_s/r)*u_phi[i][qp]*jac;

            if( this->_flow_vars.dim() > 1 )
              (*Fv)(i) += ( - tau_C*RC_s*u_gradphi[i][qp](1)
                            - tau_M*RM_s(1)*rho*U*u_gradphi[i][qp] )*jac;

            if( this->_flow_vars.dim() == 3 )
              (*Fw)(i) += ( - tau_C*RC_s*u_gradphi[i][qp](2)
                            - tau_M*RM_s(2)*rho*U*u_gradphi[i][qp] )*jac;
          }

        // Energy SUPG terms
        for (unsigned int i=0; i != n_T_dofs; i++)
          FT(i) -= rho*cp*tau_E*RE_s*U*T_gradphi[i][qp]*jac;

        // Species SUPG terms
        for(unsigned int s=0; s < this->n_species(); s++)
          {
            libMesh::DenseSubVector<libMesh::Number> &Fs =
              context.get_elem_residual(this->_species_vars.species(s));

            libMesh::Real tau_s = this->_stab_helper.compute_tau_species( context, qp, g, G, rho, U, D[s], this->_is_steady );

            for (unsigned int i=0; i != n_s_dofs; i++)
              Fs(i) -= rho*tau_s*Rs_s[s]*U*s_gradphi[i][qp]*jac;
          }

        if(compute_jacobian)
          libmesh_not_implemented();
      }
  }

  template<typename Mixture, typename Evaluator>
  void ReactingLowMachNavierStokesSPGSMStabilization<Mixture,Evaluator>::mass_residual
  ( bool compute_jacobian, AssemblyContext & context )
  {
    const unsigned int n_p_dofs = context.get_dof_indices(this->_press_var.p()).size();
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();
    const VariableIndex s0_var = this->_species_vars.species(0);
    const unsigned int n_s_dofs = context.get_dof_indices(s0_var).size();

    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
      context.get_element_fe(this->_press_var.p())->get_dphi();

    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_flow_vars.u())->get_dphi();

    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(this->_temp_vars.T())->get_dphi();

    const std::vector<std::vector<libMesh::Gradient> >& s_gradphi = context.get_element_fe(s0_var)->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(this->_press_var.p()); // R_{p}
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u()); // R_{u}

    libMesh::DenseSubVector<libMesh::Number>* Fv = NULL;
    if( this->_flow_vars.dim() > 1 )
      Fv  = &context.get_elem_residual(this->_flow_vars.v()); // R_{v}

    libMesh::DenseSubVector<libMesh::Number>* Fw = NULL;
    if( this->_flow_vars.dim() == 3 )
      Fw  = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}

    libMesh::DenseSubVector<libMesh::Number> &FT = context.get_elem_residual(this->_temp_vars.T());

    libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    const std::vector<libMesh::Point>& u_qpoint =
      context.get_element_fe(this->_flow_vars.u())->get_xyz();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::Real T = context.interior_value( this->_temp_vars.T(), qp );

        libMesh::Number u = context.fixed_interior_value(this->_flow_vars.u(), qp);

        libMesh::NumberVectorValue U(u);
        if (this->_flow_vars.dim() > 1)
          U(1) = context.fixed_interior_value(this->_flow_vars.v(), qp);
        if (this->_flow_vars.dim() == 3)
          U(2) = context.fixed_interior_value(this->_flow_vars.w(), qp);

        std::vector<libMesh::Real> ws(this->n_species());
        for(unsigned int s=0; s < this->_n_species; s++ )
          ws[s] = context.fixed_interior_value(this->_species_vars.species(s), qp);

        Evaluator gas_evaluator( *(this->_gas_mixture) );
        const libMesh::Real R_mix = gas_evaluator.R_mix(ws);
        const libMesh::Real p0 = this->get_p0_steady(context,qp);
        libMesh::Real rho = this->rho(T, p0, R_mix);

        const libMesh::Real cp = gas_evaluator.cp(T,p0,ws);

        std::vector<libMesh::Real> D( this->n_species() );
        libMesh::Real mu, k;

        gas_evaluator.mu_and_k_and_D( T, rho, cp, ws, mu, k, D );

        libMesh::Real tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, rho, U, mu, false );
        libMesh::Real tau_C = this->_stab_helper.compute_tau_continuity( tau_M, g );
        libMesh::Real tau_E = this->_stab_helper.compute_tau_energy( context, qp, g, G, rho, U, k, cp, false );

        // Strong residuals
        libMesh::Real RC_t;
        libMesh::RealGradient RM_t;
        libMesh::Real RE_t;
        std::vector<libMesh::Real> Rs_t(this->n_species());

        this->compute_res_transient( context, qp, RC_t, RM_t, RE_t, Rs_t );

        libMesh::Real jac = JxW[qp];
        const libMesh::Number r = u_qpoint[qp](0);

        if( this->_is_axisymmetric )
          jac *= r;

        for (unsigned int i=0; i != n_p_dofs; i++)
          Fp(i) -= tau_M*RM_t*p_dphi[i][qp]*jac;

        for(unsigned int s=0; s < this->n_species(); s++)
          {
            libMesh::DenseSubVector<libMesh::Number> &Fs =
              context.get_elem_residual(this->_species_vars.species(s));

            libMesh::Real tau_s = this->_stab_helper.compute_tau_species( context, qp, g, G, rho, U, D[s], false );

            for (unsigned int i=0; i != n_s_dofs; i++)
              Fs(i) -= rho*tau_s*Rs_t[s]*U*s_gradphi[i][qp]*jac;
          }

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) -= ( tau_C*RC_t*u_gradphi[i][qp](0)
                       + tau_M*RM_t(0)*rho*U*u_gradphi[i][qp] )*jac;

            if( this->_is_axisymmetric )
              Fu(i) += (-tau_C*RC_t/r)*u_phi[i][qp]*jac;

            if( this->_flow_vars.dim() > 1 )
              (*Fv)(i) -= ( tau_C*RC_t*u_gradphi[i][qp](1)
                            + tau_M*RM_t(1)*rho*U*u_gradphi[i][qp] )*jac;

            if( this->_flow_vars.dim() == 3 )
              (*Fw)(i) -= ( tau_C*RC_t*u_gradphi[i][qp](2)
                            + tau_M*RM_t(2)*rho*U*u_gradphi[i][qp] )*jac;
          }

        for (unsigned int i=0; i != n_T_dofs; i++)
          FT(i) -= rho*cp*tau_E*RE_t*U*T_gradphi[i][qp]*jac;

        if(compute_jacobian)
          libmesh_not_implemented();
      }

    return;
  }
} // end namespace GRINS
