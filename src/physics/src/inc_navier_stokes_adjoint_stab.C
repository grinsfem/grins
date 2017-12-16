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
#include "grins/inc_navier_stokes_adjoint_stab.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/inc_nav_stokes_macro.h"

//libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{

  template<class Mu>
  IncompressibleNavierStokesAdjointStabilization<Mu>::IncompressibleNavierStokesAdjointStabilization( const std::string& physics_name,
                                                                                                      const GetPot& input )
    : IncompressibleNavierStokesStabilizationBase<Mu>(physics_name,input)
  {}

  template<class Mu>
  void IncompressibleNavierStokesAdjointStabilization<Mu>::element_time_derivative
  ( bool compute_jacobian,
    AssemblyContext & context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.get_dof_indices(this->_press_var.p()).size();
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    const std::vector<std::vector<libMesh::RealGradient> >& p_gradphi =
      context.get_element_fe(this->_press_var.p())->get_dphi();

    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_flow_vars.u())->get_dphi();

    const std::vector<std::vector<libMesh::RealTensor> >& u_hessphi =
      context.get_element_fe(this->_flow_vars.u())->get_d2phi();

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u()); // R_{p}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v()); // R_{p}
    libMesh::DenseSubMatrix<libMesh::Number> &Kup =
      context.get_elem_jacobian(this->_flow_vars.u(), this->_press_var.p()); // J_{up}
    libMesh::DenseSubMatrix<libMesh::Number> &Kuu =
      context.get_elem_jacobian(this->_flow_vars.u(), this->_flow_vars.u()); // J_{uu}
    libMesh::DenseSubMatrix<libMesh::Number> &Kuv =
      context.get_elem_jacobian(this->_flow_vars.u(), this->_flow_vars.v()); // J_{uv}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvp =
      context.get_elem_jacobian(this->_flow_vars.v(), this->_press_var.p()); // J_{vp}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvu =
      context.get_elem_jacobian(this->_flow_vars.v(), this->_flow_vars.u()); // J_{vu}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvv =
      context.get_elem_jacobian(this->_flow_vars.v(), this->_flow_vars.v()); // J_{vv}

    libMesh::DenseSubVector<libMesh::Number> *Fw = NULL;
    libMesh::DenseSubMatrix<libMesh::Number> *Kuw = NULL;
    libMesh::DenseSubMatrix<libMesh::Number> *Kvw = NULL;
    libMesh::DenseSubMatrix<libMesh::Number> *Kwp = NULL;
    libMesh::DenseSubMatrix<libMesh::Number> *Kwu = NULL;
    libMesh::DenseSubMatrix<libMesh::Number> *Kwv = NULL;
    libMesh::DenseSubMatrix<libMesh::Number> *Kww = NULL;


    if(this->_flow_vars.dim() == 3)
      {
        Fw = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}
        Kuw = &context.get_elem_jacobian
          (this->_flow_vars.u(), this->_flow_vars.w()); // J_{uw}
        Kvw = &context.get_elem_jacobian
          (this->_flow_vars.v(), this->_flow_vars.w()); // J_{vw}
        Kwp = &context.get_elem_jacobian
          (this->_flow_vars.w(), this->_press_var.p()); // J_{wp}
        Kwu = &context.get_elem_jacobian
          (this->_flow_vars.w(), this->_flow_vars.u()); // J_{wu}
        Kwv = &context.get_elem_jacobian
          (this->_flow_vars.w(), this->_flow_vars.v()); // J_{wv}
        Kww = &context.get_elem_jacobian
          (this->_flow_vars.w(), this->_flow_vars.w()); // J_{ww}
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u());

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::RealGradient U( context.interior_value( this->_flow_vars.u(), qp ),
                                 context.interior_value( this->_flow_vars.v(), qp ) );
        if( this->_flow_vars.dim() == 3 )
          {
            U(2) = context.interior_value( this->_flow_vars.w(), qp );
          }

        /*
          libMesh::Gradient grad_u, grad_v, grad_w;
          grad_u = context.interior_gradient(this->_flow_vars.u(), qp);
          grad_v = context.interior_gradient(this->_flow_vars.v(), qp);
          if (_dim == 3)
          grad_w = context.interior_gradient(this->_flow_vars.w(), qp);
        */

        libMesh::Real tau_M, tau_C;
        libMesh::Real d_tau_M_d_rho, d_tau_C_d_rho;
        libMesh::Gradient d_tau_M_dU, d_tau_C_dU;
        libMesh::Gradient RM_s, d_RM_s_uvw_dgraduvw;
        libMesh::Real RC;
        libMesh::Tensor d_RC_dgradU,
          d_RM_s_dgradp, d_RM_s_dU, d_RM_s_uvw_dhessuvw;

        // Compute the viscosity at this qp
        libMesh::Real _mu_qp = this->_mu(context, qp);

        if (compute_jacobian)
          {
            this->_stab_helper.compute_tau_momentum_and_derivs
              ( context, qp, g, G, this->_rho, U, _mu_qp,
                tau_M, d_tau_M_d_rho, d_tau_M_dU,
                this->_is_steady );
            this->_stab_helper.compute_tau_continuity_and_derivs
              ( tau_M, d_tau_M_d_rho, d_tau_M_dU,
                g,
                tau_C, d_tau_C_d_rho, d_tau_C_dU );
            this->_stab_helper.compute_res_momentum_steady_and_derivs
              ( context, qp, this->_rho, _mu_qp,
                RM_s, d_RM_s_dgradp, d_RM_s_dU, d_RM_s_uvw_dgraduvw,
                d_RM_s_uvw_dhessuvw);
            this->_stab_helper.compute_res_continuity_and_derivs
              ( context, qp, RC, d_RC_dgradU );
          }
        else
          {
            tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, this->_rho, U, _mu_qp, this->_is_steady );
            tau_C = this->_stab_helper.compute_tau_continuity( tau_M, g );
            RM_s = this->_stab_helper.compute_res_momentum_steady( context, qp, this->_rho, _mu_qp );
            RC = this->_stab_helper.compute_res_continuity( context, qp );
          }

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            libMesh::Real test_func = this->_rho*U*u_gradphi[i][qp] +
              _mu_qp*( u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) + u_hessphi[i][qp](2,2) );
            libMesh::Gradient d_test_func_dU = this->_rho*u_gradphi[i][qp];

            //libMesh::RealGradient zeroth_order_term = - this->_rho*u_phi[i][qp]*(grad_u + grad_v + grad_w);

            Fu(i) += ( -tau_M*RM_s(0)*test_func - tau_C*RC*u_gradphi[i][qp](0) )*JxW[qp];

            Fv(i) += ( -tau_M*RM_s(1)*test_func - tau_C*RC*u_gradphi[i][qp](1) )*JxW[qp];

            if(this->_flow_vars.dim() == 3)
              {
                (*Fw)(i) += ( -tau_M*RM_s(2)*test_func - tau_C*RC*u_gradphi[i][qp](2) )*JxW[qp];
              }

            if (compute_jacobian)
              {
                const libMesh::Real fixed_deriv =
                  context.get_fixed_solution_derivative();
                for (unsigned int j=0; j != n_p_dofs; j++)
                  {
                    Kup(i,j) += ( -tau_M*
                                  ( d_RM_s_dgradp(0,0) *
                                    p_gradphi[j][qp](0) +
                                    d_RM_s_dgradp(0,1) *
                                    p_gradphi[j][qp](1) +
                                    d_RM_s_dgradp(0,2) *
                                    p_gradphi[j][qp](2)
                                    )*test_func)*fixed_deriv*JxW[qp];
                    Kvp(i,j) += ( -tau_M*
                                  ( d_RM_s_dgradp(1,0) *
                                    p_gradphi[j][qp](0) +
                                    d_RM_s_dgradp(1,1) *
                                    p_gradphi[j][qp](1) +
                                    d_RM_s_dgradp(1,2) *
                                    p_gradphi[j][qp](2)
                                    )*test_func)*fixed_deriv*JxW[qp];
                  }
                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    Kuu(i,j) += ( -d_tau_M_dU(0)*RM_s(0)*test_func
                                  -tau_M*d_RM_s_dU(0,0)*test_func
                                  - d_tau_C_dU(0)*RC*u_gradphi[i][qp](0)
                                  )*fixed_deriv*u_phi[j][qp]*JxW[qp];
                    Kuu(i,j) += ( -tau_M*RM_s(0)*d_test_func_dU(0)
                                  )*u_phi[j][qp]*JxW[qp];
                    Kuu(i,j) += ( - tau_C*
                                  (
                                   d_RC_dgradU(0,0)*u_gradphi[j][qp](0) +
                                   d_RC_dgradU(0,1)*u_gradphi[j][qp](1) +
                                   d_RC_dgradU(0,2)*u_gradphi[j][qp](2)
                                   )
                                  *fixed_deriv*u_gradphi[i][qp](0)
                                  )*JxW[qp];
                    Kuu(i,j) += ( -tau_M*test_func*(d_RM_s_uvw_dgraduvw*u_gradphi[j][qp])
                                  )*fixed_deriv*JxW[qp];
                    Kuu(i,j) += ( -tau_M*test_func*(d_RM_s_uvw_dhessuvw.contract(u_hessphi[j][qp]))
                                  )*fixed_deriv*JxW[qp];
                    Kuv(i,j) += ( -d_tau_M_dU(1)*RM_s(0)*test_func
                                  -tau_M*d_RM_s_dU(0,1)*test_func
                                  )*fixed_deriv*u_phi[j][qp]*JxW[qp];
                    Kuv(i,j) += ( -tau_M*RM_s(0)*d_test_func_dU(1)
                                  )*u_phi[j][qp]*JxW[qp];
                    Kuv(i,j) += ( - tau_C*
                                  (
                                   d_RC_dgradU(1,0)*u_gradphi[j][qp](0) +
                                   d_RC_dgradU(1,1)*u_gradphi[j][qp](1) +
                                   d_RC_dgradU(1,2)*u_gradphi[j][qp](2)
                                   )
                                  *fixed_deriv*u_gradphi[i][qp](0)
                                  )*JxW[qp];
                    Kvu(i,j) += ( -d_tau_M_dU(0)*RM_s(1)*test_func
                                  -tau_M*d_RM_s_dU(1,0)*test_func
                                  )*fixed_deriv*u_phi[j][qp]*JxW[qp];
                    Kvu(i,j) += ( -tau_M*RM_s(1)*d_test_func_dU(0)
                                  )*u_phi[j][qp]*JxW[qp];
                    Kvu(i,j) += ( - tau_C*
                                  (
                                   d_RC_dgradU(0,0)*u_gradphi[j][qp](0) +
                                   d_RC_dgradU(0,1)*u_gradphi[j][qp](1) +
                                   d_RC_dgradU(0,2)*u_gradphi[j][qp](2)
                                   )
                                  *fixed_deriv*u_gradphi[i][qp](1)
                                  )*JxW[qp];
                    Kvv(i,j) += ( -d_tau_M_dU(1)*RM_s(1)*test_func
                                  -tau_M*d_RM_s_dU(1,1)*test_func
                                  )*fixed_deriv*u_phi[j][qp]*JxW[qp];
                    Kvv(i,j) += ( -tau_M*RM_s(1)*d_test_func_dU(1)
                                  )*u_phi[j][qp]*JxW[qp];
                    Kvv(i,j) += ( - tau_C*
                                  (
                                   d_RC_dgradU(1,0)*u_gradphi[j][qp](0) +
                                   d_RC_dgradU(1,1)*u_gradphi[j][qp](1) +
                                   d_RC_dgradU(1,2)*u_gradphi[j][qp](2)
                                   )
                                  *fixed_deriv*u_gradphi[i][qp](1)
                                  )*JxW[qp];
                    Kvv(i,j) += ( -tau_M*test_func*(d_RM_s_uvw_dgraduvw*u_gradphi[j][qp])
                                  )*fixed_deriv*JxW[qp];
                    Kvv(i,j) += ( -tau_M*test_func*(d_RM_s_uvw_dhessuvw.contract(u_hessphi[j][qp]))
                                  )*fixed_deriv*JxW[qp];
                  }
                if(this->_flow_vars.dim() == 3)
                  {
                    for (unsigned int j=0; j != n_p_dofs; j++)
                      {
                        (*Kwp)(i,j) += ( -tau_M*
                                         ( d_RM_s_dgradp(2,0) *
                                           p_gradphi[j][qp](0) +
                                           d_RM_s_dgradp(2,1) *
                                           p_gradphi[j][qp](1) +
                                           d_RM_s_dgradp(2,2) *
                                           p_gradphi[j][qp](2)
                                           )*test_func)*fixed_deriv*JxW[qp];
                      }
                    for (unsigned int j=0; j != n_u_dofs; j++)
                      {
                        (*Kuw)(i,j) += ( -d_tau_M_dU(2)*RM_s(0)*test_func
                                         -tau_M*d_RM_s_dU(0,2)*test_func
                                         )*fixed_deriv*u_phi[j][qp]*JxW[qp];
                        (*Kuw)(i,j) += ( -tau_M*RM_s(0)*d_test_func_dU(2)
                                         )*u_phi[j][qp]*JxW[qp];
                        (*Kuw)(i,j) += ( - tau_C*
                                         (
                                          d_RC_dgradU(2,0)*u_gradphi[j][qp](0) +
                                          d_RC_dgradU(2,1)*u_gradphi[j][qp](1) +
                                          d_RC_dgradU(2,2)*u_gradphi[j][qp](2)
                                          )
                                         *fixed_deriv* u_gradphi[i][qp](0)
                                         )*JxW[qp];
                        (*Kvw)(i,j) += ( -d_tau_M_dU(2)*RM_s(1)*test_func
                                         -tau_M*d_RM_s_dU(1,2)*test_func
                                         )*fixed_deriv*u_phi[j][qp]*JxW[qp];
                        (*Kvw)(i,j) += ( -tau_M*RM_s(1)*d_test_func_dU(2)
                                         )*u_phi[j][qp]*JxW[qp];
                        (*Kvw)(i,j) += ( - tau_C*
                                         (
                                          d_RC_dgradU(2,0)*u_gradphi[j][qp](0) +
                                          d_RC_dgradU(2,1)*u_gradphi[j][qp](1) +
                                          d_RC_dgradU(2,2)*u_gradphi[j][qp](2)
                                          )
                                         *fixed_deriv* u_gradphi[i][qp](1)
                                         )*JxW[qp];
                        (*Kwu)(i,j) += ( -d_tau_M_dU(0)*RM_s(2)*test_func
                                         -tau_M*d_RM_s_dU(2,0)*test_func
                                         )*fixed_deriv*u_phi[j][qp]*JxW[qp];
                        (*Kwu)(i,j) += ( -tau_M*RM_s(2)*d_test_func_dU(0)
                                         )*u_phi[j][qp]*JxW[qp];
                        (*Kwu)(i,j) += ( - tau_C*
                                         (
                                          d_RC_dgradU(0,0)*u_gradphi[j][qp](0) +
                                          d_RC_dgradU(0,1)*u_gradphi[j][qp](1) +
                                          d_RC_dgradU(0,2)*u_gradphi[j][qp](2)
                                          )
                                         *fixed_deriv* u_gradphi[i][qp](2)
                                         )*JxW[qp];
                        (*Kwv)(i,j) += ( -d_tau_M_dU(1)*RM_s(2)*test_func
                                         -tau_M*d_RM_s_dU(2,1)*test_func
                                         )*fixed_deriv*u_phi[j][qp]*JxW[qp];
                        (*Kwv)(i,j) += ( -tau_M*RM_s(2)*d_test_func_dU(1)
                                         )*u_phi[j][qp]*JxW[qp];
                        (*Kwv)(i,j) += ( - tau_C*
                                         (
                                          d_RC_dgradU(1,0)*u_gradphi[j][qp](0) +
                                          d_RC_dgradU(1,1)*u_gradphi[j][qp](1) +
                                          d_RC_dgradU(1,2)*u_gradphi[j][qp](2)
                                          )
                                         *fixed_deriv* u_gradphi[i][qp](2)
                                         )*JxW[qp];
                        (*Kww)(i,j) += ( -d_tau_M_dU(2)*RM_s(2)*test_func
                                         -tau_M*d_RM_s_dU(2,2)*test_func
                                         )*fixed_deriv*u_phi[j][qp]*JxW[qp];
                        (*Kww)(i,j) += ( -tau_M*RM_s(2)*d_test_func_dU(2)
                                         )*u_phi[j][qp]*JxW[qp];
                        (*Kww)(i,j) += ( - tau_C*
                                         (
                                          d_RC_dgradU(2,0)*u_gradphi[j][qp](0) +
                                          d_RC_dgradU(2,1)*u_gradphi[j][qp](1) +
                                          d_RC_dgradU(2,2)*u_gradphi[j][qp](2)
                                          )
                                         *fixed_deriv* u_gradphi[i][qp](2)
                                         )*JxW[qp];
                        (*Kww)(i,j) += ( -tau_M*test_func*(d_RM_s_uvw_dgraduvw*u_gradphi[j][qp])
                                         )*fixed_deriv*JxW[qp];
                        (*Kww)(i,j) += ( -tau_M*test_func*(d_RM_s_uvw_dhessuvw.contract(u_hessphi[j][qp]))
                                         )*fixed_deriv*JxW[qp];

                      }
                  }
              }
          }
      }
  }

  template<class Mu>
  void IncompressibleNavierStokesAdjointStabilization<Mu>::element_constraint
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.get_dof_indices(this->_press_var.p()).size();
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
      context.get_element_fe(this->_press_var.p())->get_dphi();

    // The velocity shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> >& u_dphi =
      context.get_element_fe(this->_flow_vars.u())->get_dphi();

    const std::vector<std::vector<libMesh::RealTensor> >& u_d2phi =
      context.get_element_fe(this->_flow_vars.u())->get_d2phi();

    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(this->_press_var.p()); // R_{p}

    libMesh::DenseSubMatrix<libMesh::Number> &Kpp =
      context.get_elem_jacobian(this->_press_var.p(), this->_press_var.p()); // J_{pp}
    libMesh::DenseSubMatrix<libMesh::Number> &Kpu =
      context.get_elem_jacobian(this->_press_var.p(), this->_flow_vars.u()); // J_{pu}
    libMesh::DenseSubMatrix<libMesh::Number> &Kpv =
      context.get_elem_jacobian(this->_press_var.p(), this->_flow_vars.v()); // J_{pv}
    libMesh::DenseSubMatrix<libMesh::Number> *Kpw = NULL;


    if(this->_flow_vars.dim() == 3)
      {
        Kpw = &context.get_elem_jacobian
          (this->_press_var.p(), this->_flow_vars.w()); // J_{pw}
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u());

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::RealGradient U( context.interior_value( this->_flow_vars.u(), qp ),
                                 context.interior_value( this->_flow_vars.v(), qp ) );
        if( this->_flow_vars.dim() == 3 )
          {
            U(2) = context.interior_value( this->_flow_vars.w(), qp );
          }

        libMesh::Real tau_M;
        libMesh::Real d_tau_M_d_rho;
        libMesh::Gradient d_tau_M_dU;
        libMesh::Gradient RM_s, d_RM_s_uvw_dgraduvw;
        libMesh::Tensor d_RM_s_dgradp, d_RM_s_dU, d_RM_s_uvw_dhessuvw;

        // Compute the viscosity at this qp
        libMesh::Real _mu_qp = this->_mu(context, qp);

        if (compute_jacobian)
          {
            this->_stab_helper.compute_tau_momentum_and_derivs
              ( context, qp, g, G, this->_rho, U, _mu_qp,
                tau_M, d_tau_M_d_rho, d_tau_M_dU,
                this->_is_steady );
            this->_stab_helper.compute_res_momentum_steady_and_derivs
              ( context, qp, this->_rho, _mu_qp,
                RM_s, d_RM_s_dgradp, d_RM_s_dU, d_RM_s_uvw_dgraduvw,
                d_RM_s_uvw_dhessuvw);
          }
        else
          {
            tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, this->_rho, U, _mu_qp, this->_is_steady );
            RM_s = this->_stab_helper.compute_res_momentum_steady( context, qp, this->_rho, _mu_qp );
          }

        // Now a loop over the pressure degrees of freedom.  This
        // computes the contributions of the continuity equation.
        for (unsigned int i=0; i != n_p_dofs; i++)
          {
            Fp(i) += tau_M*(RM_s*p_dphi[i][qp])*JxW[qp];

            if (compute_jacobian)
              {
                const libMesh::Real fixed_deriv =
                  context.get_fixed_solution_derivative();

                const libMesh::TypeVector<libMesh::Number>
                  p_dphiiT_d_RM_s_dgradp =
                  d_RM_s_dgradp.transpose() * p_dphi[i][qp];

                for (unsigned int j=0; j != n_p_dofs; j++)
                  {
                    Kpp(i,j) += tau_M*(p_dphiiT_d_RM_s_dgradp*p_dphi[j][qp])*JxW[qp];
                  }

                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    Kpu(i,j) += (
                                 d_tau_M_dU(0)*u_phi[j][qp]*(RM_s*p_dphi[i][qp])
                                 + tau_M*
                                 (
                                  d_RM_s_dU(0,0)*u_phi[j][qp]*p_dphi[i][qp](0) +
                                  d_RM_s_dU(1,0)*u_phi[j][qp]*p_dphi[i][qp](1) +
                                  d_RM_s_dU(2,0)*u_phi[j][qp]*p_dphi[i][qp](2) +
                                  d_RM_s_uvw_dgraduvw*u_dphi[j][qp]*p_dphi[i][qp](0) +
                                  d_RM_s_uvw_dhessuvw.contract(u_d2phi[j][qp])*p_dphi[i][qp](0)
                                  )*fixed_deriv
                                 )*JxW[qp];
                    Kpv(i,j) += (
                                 d_tau_M_dU(1)*u_phi[j][qp]*(RM_s*p_dphi[i][qp])
                                 + tau_M*
                                 (
                                  d_RM_s_dU(0,1)*u_phi[j][qp]*p_dphi[i][qp](0) +
                                  d_RM_s_dU(1,1)*u_phi[j][qp]*p_dphi[i][qp](1) +
                                  d_RM_s_dU(2,1)*u_phi[j][qp]*p_dphi[i][qp](2) +
                                  d_RM_s_uvw_dgraduvw*u_dphi[j][qp]*p_dphi[i][qp](1) +
                                  d_RM_s_uvw_dhessuvw.contract(u_d2phi[j][qp])*p_dphi[i][qp](0)
                                  )*fixed_deriv
                                 )*JxW[qp];
                  }

                if(this->_flow_vars.dim() == 3)
                  {
                    for (unsigned int j=0; j != n_u_dofs; j++)
                      {
                        (*Kpw)(i,j) += (
                                        d_tau_M_dU(2)*u_phi[j][qp]*(RM_s*p_dphi[i][qp])
                                        + tau_M*
                                        (
                                         d_RM_s_dU(0,2)*u_phi[j][qp]*p_dphi[i][qp](0) +
                                         d_RM_s_dU(1,2)*u_phi[j][qp]*p_dphi[i][qp](1) +
                                         d_RM_s_dU(2,2)*u_phi[j][qp]*p_dphi[i][qp](2) +
                                         d_RM_s_uvw_dgraduvw*u_dphi[j][qp]*p_dphi[i][qp](2) +
                                         d_RM_s_uvw_dhessuvw.contract(u_d2phi[j][qp])*p_dphi[i][qp](0)
                                         )*fixed_deriv
                                        )*JxW[qp];
                      }
                  }
              }
          }
      }
  }

  template<class Mu>
  void IncompressibleNavierStokesAdjointStabilization<Mu>::mass_residual
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.get_dof_indices(this->_press_var.p()).size();
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();

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

    const std::vector<std::vector<libMesh::RealTensor> >& u_hessphi =
      context.get_element_fe(this->_flow_vars.u())->get_d2phi();

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u()); // R_{p}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v()); // R_{p}
    libMesh::DenseSubVector<libMesh::Number> *Fw = NULL;

    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(this->_press_var.p()); // R_{p}

    libMesh::DenseSubMatrix<libMesh::Number> &Kuu =
      context.get_elem_jacobian(this->_flow_vars.u(), this->_flow_vars.u()); // J_{uu}
    libMesh::DenseSubMatrix<libMesh::Number> &Kuv =
      context.get_elem_jacobian(this->_flow_vars.u(), this->_flow_vars.v()); // J_{uv}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvu =
      context.get_elem_jacobian(this->_flow_vars.v(), this->_flow_vars.u()); // J_{vu}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvv =
      context.get_elem_jacobian(this->_flow_vars.v(), this->_flow_vars.v()); // J_{vv}

    libMesh::DenseSubMatrix<libMesh::Number> *Kuw = NULL;
    libMesh::DenseSubMatrix<libMesh::Number> *Kvw = NULL;
    libMesh::DenseSubMatrix<libMesh::Number> *Kwu = NULL;
    libMesh::DenseSubMatrix<libMesh::Number> *Kwv = NULL;
    libMesh::DenseSubMatrix<libMesh::Number> *Kww = NULL;

    libMesh::DenseSubMatrix<libMesh::Number> &Kpu =
      context.get_elem_jacobian(this->_press_var.p(), this->_flow_vars.u()); // J_{pu}
    libMesh::DenseSubMatrix<libMesh::Number> &Kpv =
      context.get_elem_jacobian(this->_press_var.p(), this->_flow_vars.v()); // J_{pv}
    libMesh::DenseSubMatrix<libMesh::Number> *Kpw = NULL;


    if(this->_flow_vars.dim() == 3)
      {
        Fw = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}

        Kuw = &context.get_elem_jacobian
          (this->_flow_vars.u(), this->_flow_vars.w()); // J_{uw}
        Kvw = &context.get_elem_jacobian
          (this->_flow_vars.v(), this->_flow_vars.w()); // J_{vw}
        Kwu = &context.get_elem_jacobian
          (this->_flow_vars.w(), this->_flow_vars.u()); // J_{wu}
        Kwv = &context.get_elem_jacobian
          (this->_flow_vars.w(), this->_flow_vars.v()); // J_{wv}
        Kww = &context.get_elem_jacobian
          (this->_flow_vars.w(), this->_flow_vars.w()); // J_{ww}

        Kpw = &context.get_elem_jacobian
          (this->_press_var.p(), this->_flow_vars.w()); // J_{pw}
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u());

        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::RealGradient U( context.fixed_interior_value( this->_flow_vars.u(), qp ),
                                 context.fixed_interior_value( this->_flow_vars.v(), qp ) );
        if( this->_flow_vars.dim() == 3 )
          {
            U(2) = context.fixed_interior_value( this->_flow_vars.w(), qp );
          }

        /*
          libMesh::Gradient grad_u, grad_v, grad_w;
          grad_u = context.interior_gradient(this->_flow_vars.u(), qp);
          grad_v = context.interior_gradient(this->_flow_vars.v(), qp);
          if (_dim == 3)
          grad_w = context.interior_gradient(this->_flow_vars.w(), qp);
        */

        libMesh::Real tau_M;

        libMesh::Real d_tau_M_d_rho;

        libMesh::Gradient d_tau_M_dU;

        libMesh::RealGradient RM_t;

        libMesh::Real d_RM_t_uvw_duvw;

        // Compute the viscosity at this qp
        libMesh::Real _mu_qp = this->_mu(context, qp);

        if (compute_jacobian)
          {
            this->_stab_helper.compute_tau_momentum_and_derivs
              ( context, qp, g, G, this->_rho, U, _mu_qp,
                tau_M, d_tau_M_d_rho, d_tau_M_dU,
                false );
            this->_stab_helper.compute_res_momentum_transient_and_derivs
              ( context, qp, this->_rho,
                RM_t, d_RM_t_uvw_duvw );
          }
        else
          {
            tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, this->_rho, U, _mu_qp, false );
            RM_t = this->_stab_helper.compute_res_momentum_transient( context, qp, this->_rho );
          }


        // Now a loop over the pressure degrees of freedom.  This
        // computes the contributions of the continuity equation.
        for (unsigned int i=0; i != n_p_dofs; i++)
          {
            Fp(i) -= tau_M*RM_t*p_dphi[i][qp]*JxW[qp];

            if (compute_jacobian)
              {
                const libMesh::Real fixed_deriv =
                  context.get_fixed_solution_derivative();

                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    Kpu(i,j) -= d_tau_M_dU(0)*u_phi[j][qp]*RM_t*p_dphi[i][qp]*fixed_deriv*JxW[qp];
                    Kpu(i,j) -= tau_M*d_RM_t_uvw_duvw*u_phi[j][qp]*p_dphi[i][qp](0)*fixed_deriv*JxW[qp];

                    Kpv(i,j) -= d_tau_M_dU(1)*u_phi[j][qp]*RM_t*p_dphi[i][qp]*fixed_deriv*JxW[qp];
                    Kpv(i,j) -= tau_M*d_RM_t_uvw_duvw*u_phi[j][qp]*p_dphi[i][qp](1)*fixed_deriv*JxW[qp];
                  }

                if(this->_flow_vars.dim() == 3)
                  {
                    for (unsigned int j=0; j != n_u_dofs; j++)
                      {
                        (*Kpw)(i,j) -= d_tau_M_dU(2)*u_phi[j][qp]*RM_t*p_dphi[i][qp]*fixed_deriv*JxW[qp];
                        (*Kpw)(i,j) -= tau_M*d_RM_t_uvw_duvw*u_phi[j][qp]*p_dphi[i][qp](2)*fixed_deriv*JxW[qp];
                      }
                  }
              }
          }

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            libMesh::Real test_func = this->_rho*U*u_gradphi[i][qp] +
              _mu_qp*( u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) + u_hessphi[i][qp](2,2) );
            libMesh::Gradient d_test_func_dU = this->_rho*u_gradphi[i][qp];

            //libMesh::RealGradient zeroth_order_term = - this->_rho*u_phi[i][qp]*(grad_u + grad_v + grad_w);

            Fu(i) -= tau_M*RM_t(0)*test_func*JxW[qp];

            Fv(i) -= tau_M*RM_t(1)*test_func*JxW[qp];

            if(this->_flow_vars.dim() == 3)
              {
                (*Fw)(i) -= tau_M*RM_t(2)*test_func*JxW[qp];
              }

            if (compute_jacobian)
              {
                const libMesh::Real fixed_deriv =
                  context.get_fixed_solution_derivative();

                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    Kuu(i,j) -= d_tau_M_dU(0)*u_phi[j][qp]*RM_t(0)*test_func*fixed_deriv*JxW[qp];
                    Kuu(i,j) -= tau_M*d_RM_t_uvw_duvw*u_phi[j][qp]*test_func*fixed_deriv*JxW[qp];
                    Kuu(i,j) -= tau_M*RM_t(0)*d_test_func_dU(0)*u_phi[j][qp]*fixed_deriv*JxW[qp];

                    Kuv(i,j) -= d_tau_M_dU(1)*u_phi[j][qp]*RM_t(0)*test_func*fixed_deriv*JxW[qp];
                    Kuv(i,j) -= tau_M*RM_t(0)*d_test_func_dU(1)*u_phi[j][qp]*fixed_deriv*JxW[qp];

                    Kvu(i,j) -= d_tau_M_dU(0)*u_phi[j][qp]*RM_t(1)*test_func*fixed_deriv*JxW[qp];
                    Kvu(i,j) -= tau_M*RM_t(1)*d_test_func_dU(0)*u_phi[j][qp]*fixed_deriv*JxW[qp];

                    Kvv(i,j) -= d_tau_M_dU(1)*u_phi[j][qp]*RM_t(1)*test_func*fixed_deriv*JxW[qp];
                    Kvv(i,j) -= tau_M*d_RM_t_uvw_duvw*u_phi[j][qp]*test_func*fixed_deriv*JxW[qp];
                    Kvv(i,j) -= tau_M*RM_t(1)*d_test_func_dU(1)*u_phi[j][qp]*fixed_deriv*JxW[qp];
                  }
                if(this->_flow_vars.dim() == 3)
                  {
                    for (unsigned int j=0; j != n_u_dofs; j++)
                      {
                        (*Kuw)(i,j) -= d_tau_M_dU(2)*u_phi[j][qp]*RM_t(0)*test_func*fixed_deriv*JxW[qp];
                        (*Kuw)(i,j) -= tau_M*RM_t(0)*d_test_func_dU(2)*u_phi[j][qp]*fixed_deriv*JxW[qp];

                        (*Kvw)(i,j) -= d_tau_M_dU(2)*u_phi[j][qp]*RM_t(1)*test_func*fixed_deriv*JxW[qp];
                        (*Kvw)(i,j) -= tau_M*RM_t(1)*d_test_func_dU(2)*u_phi[j][qp]*fixed_deriv*JxW[qp];

                        (*Kwu)(i,j) -= d_tau_M_dU(0)*u_phi[j][qp]*RM_t(2)*test_func*fixed_deriv*JxW[qp];
                        (*Kwu)(i,j) -= tau_M*RM_t(2)*d_test_func_dU(0)*u_phi[j][qp]*fixed_deriv*JxW[qp];

                        (*Kwv)(i,j) -= d_tau_M_dU(1)*u_phi[j][qp]*RM_t(2)*test_func*fixed_deriv*JxW[qp];
                        (*Kwv)(i,j) -= tau_M*RM_t(2)*d_test_func_dU(1)*u_phi[j][qp]*fixed_deriv*JxW[qp];

                        (*Kww)(i,j) -= d_tau_M_dU(2)*u_phi[j][qp]*RM_t(2)*test_func*fixed_deriv*JxW[qp];
                        (*Kww)(i,j) -= tau_M*d_RM_t_uvw_duvw*u_phi[j][qp]*test_func*fixed_deriv*JxW[qp];
                        (*Kww)(i,j) -= tau_M*RM_t(2)*d_test_func_dU(2)*u_phi[j][qp]*fixed_deriv*JxW[qp];
                      }
                  }
              }
          }

      }
  }

} // namespace GRINS

// Instantiate
INSTANTIATE_INC_NS_SUBCLASS(IncompressibleNavierStokesAdjointStabilization);
