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
#include "grins/low_mach_navier_stokes_braack_stab.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/constant_viscosity.h"
#include "grins/constant_specific_heat.h"
#include "grins/constant_conductivity.h"

// libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{

  template<class Mu, class SH, class TC>
  LowMachNavierStokesBraackStabilization<Mu,SH,TC>::LowMachNavierStokesBraackStabilization( const std::string& physics_name,
                                                                                            const GetPot& input )
    : LowMachNavierStokesStabilizationBase<Mu,SH,TC>(physics_name,input)
  {
    return;
  }

  template<class Mu, class SH, class TC>
  LowMachNavierStokesBraackStabilization<Mu,SH,TC>::~LowMachNavierStokesBraackStabilization()
  {
    return;
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokesBraackStabilization<Mu,SH,TC>::element_time_derivative
  ( bool compute_jacobian,
    AssemblyContext & context )
  {
    this->assemble_continuity_time_deriv( compute_jacobian, context );
    this->assemble_momentum_time_deriv( compute_jacobian, context );
    this->assemble_energy_time_deriv( compute_jacobian, context );
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokesBraackStabilization<Mu,SH,TC>::mass_residual
  ( bool compute_jacobian, AssemblyContext & context )
  {
    this->assemble_continuity_mass_residual( compute_jacobian, context );
    this->assemble_momentum_mass_residual( compute_jacobian, context );
    this->assemble_energy_mass_residual( compute_jacobian, context );
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokesBraackStabilization<Mu,SH,TC>::assemble_continuity_time_deriv( bool /*compute_jacobian*/,
                                                                                         AssemblyContext& context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.get_dof_indices(this->_press_var.p()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
      context.get_element_fe(this->_press_var.p())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(this->_press_var.p()); // R_{p}

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u());

        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::Real T = context.interior_value( this->_temp_vars.T(), qp );
        libMesh::Real rho = this->rho( T, this->get_p0_steady( context, qp ) );

        libMesh::Real mu = this->_mu(T);
        libMesh::Real k = this->_k(T);
        libMesh::Real cp = this->_cp(T);

        libMesh::RealGradient U( context.interior_value( this->_flow_vars.u(), qp ),
                                 context.interior_value( this->_flow_vars.v(), qp ) );
        if( this->_flow_vars.dim() == 3 )
          U(2) = context.interior_value( this->_flow_vars.w(), qp ); // w

        libMesh::Real tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, rho, U, mu, this->_is_steady );
        libMesh::Real tau_E = this->_stab_helper.compute_tau_energy( context, qp, g, G, rho, U, k, cp, this->_is_steady );

        libMesh::RealGradient RM_s = this->compute_res_momentum_steady( context, qp );
        libMesh::Real RE_s = this->compute_res_energy_steady( context, qp );

        // Now a loop over the pressure degrees of freedom.  This
        // computes the contributions of the continuity equation.
        for (unsigned int i=0; i != n_p_dofs; i++)
          {
            Fp(i) += ( tau_M*RM_s*p_dphi[i][qp]
                       + tau_E*RE_s*(U*p_dphi[i][qp])/T )*JxW[qp];
          }

      }

    return;
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokesBraackStabilization<Mu,SH,TC>::assemble_momentum_time_deriv( bool /*compute_jacobian*/,
                                                                                       AssemblyContext& context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();

    // Check number of dofs is same for _flow_vars.u(), v_var and w_var.
    libmesh_assert (n_u_dofs == context.get_dof_indices(this->_flow_vars.v()).size());
    if (this->_flow_vars.dim() == 3)
      libmesh_assert (n_u_dofs == context.get_dof_indices(this->_flow_vars.w()).size());

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // The velocity shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_flow_vars.u())->get_dphi();

    const std::vector<std::vector<libMesh::RealTensor> >& u_hessphi =
      context.get_element_fe(this->_flow_vars.u())->get_d2phi();

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u()); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v()); // R_{v}
    libMesh::DenseSubVector<libMesh::Real>* Fw = NULL;

    if( this->_flow_vars.dim() == 3 )
      {
        Fw  = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Real T = context.interior_value( this->_temp_vars.T(), qp );
        libMesh::Real rho = this->rho( T, this->get_p0_steady( context, qp ) );

        libMesh::Real mu = this->_mu(T);

        libMesh::RealGradient U( context.interior_value(this->_flow_vars.u(), qp),
                                 context.interior_value(this->_flow_vars.v(), qp) );

        libMesh::RealGradient grad_u = context.interior_gradient(this->_flow_vars.u(), qp);
        libMesh::RealGradient grad_v = context.interior_gradient(this->_flow_vars.v(), qp);
        libMesh::RealGradient grad_w;

        if( this->_flow_vars.dim() == 3 )
          {
            U(2) = context.interior_value(this->_flow_vars.w(), qp);
            grad_w = context.interior_gradient(this->_flow_vars.w(), qp);
          }

        libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u());

        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::Real tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, rho, U, mu, this->_is_steady );
        libMesh::Real tau_C = this->_stab_helper.compute_tau_continuity( tau_M, g );

        libMesh::Real RC_s = this->compute_res_continuity_steady( context, qp );
        libMesh::RealGradient RM_s = this->compute_res_momentum_steady( context, qp );

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) += ( tau_C*RC_s*u_gradphi[i][qp](0)
                       + tau_M*RM_s(0)*rho*U*u_gradphi[i][qp]
                       + mu*tau_M*RM_s(0)*(u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1)
                                           + u_hessphi[i][qp](0,0) + u_hessphi[i][qp](0,1)
                                           - 2.0/3.0*(u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,0)) )
                       )*JxW[qp];

            Fv(i) += ( tau_C*RC_s*u_gradphi[i][qp](1)
                       + tau_M*RM_s(1)*rho*U*u_gradphi[i][qp]
                       + mu*tau_M*RM_s(1)*(u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1)
                                           + u_hessphi[i][qp](1,0) + u_hessphi[i][qp](1,1)
                                           - 2.0/3.0*(u_hessphi[i][qp](0,1) + u_hessphi[i][qp](1,1)) )
                       )*JxW[qp];

            if( this->_flow_vars.dim() == 3 )
              {
                Fu(i) += mu*tau_M*RM_s(0)*(u_hessphi[i][qp](2,2) + u_hessphi[i][qp](0,2)
                                           - 2.0/3.0*u_hessphi[i][qp](2,0))*JxW[qp];

                Fv(i) += mu*tau_M*RM_s(1)*(u_hessphi[i][qp](2,2) + u_hessphi[i][qp](1,2)
                                           - 2.0/3.0*u_hessphi[i][qp](2,1))*JxW[qp];

                (*Fw)(i) += ( tau_C*RC_s*u_gradphi[i][qp](2)
                              + tau_M*RM_s(2)*rho*U*u_gradphi[i][qp]
                              + mu*tau_M*RM_s(2)*(u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) + u_hessphi[i][qp](2,2)
                                                  + u_hessphi[i][qp](2,0) + u_hessphi[i][qp](2,1) + u_hessphi[i][qp](2,2)
                                                  - 2.0/3.0*(u_hessphi[i][qp](0,2) + u_hessphi[i][qp](1,2)
                                                             + u_hessphi[i][qp](2,2))
                                                  )
                              )*JxW[qp];
              }
          }

      }
    return;
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokesBraackStabilization<Mu,SH,TC>::assemble_energy_time_deriv( bool /*compute_jacobian*/,
                                                                                     AssemblyContext& context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_temp_vars.T())->get_JxW();

    // The temperature shape functions gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(this->_temp_vars.T())->get_dphi();

    const std::vector<std::vector<libMesh::RealTensor> >& T_hessphi =
      context.get_element_fe(this->_temp_vars.T())->get_d2phi();

    libMesh::DenseSubVector<libMesh::Number> &FT = context.get_elem_residual(this->_temp_vars.T()); // R_{T}

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Number u, v;
        u = context.interior_value(this->_flow_vars.u(), qp);
        v = context.interior_value(this->_flow_vars.v(), qp);

        libMesh::Gradient grad_T = context.interior_gradient(this->_temp_vars.T(), qp);

        libMesh::NumberVectorValue U(u,v);
        if (this->_flow_vars.dim() == 3)
          U(2) = context.interior_value(this->_flow_vars.w(), qp);

        libMesh::Real T = context.interior_value( this->_temp_vars.T(), qp );
        libMesh::Real rho = this->rho( T, this->get_p0_steady( context, qp ) );

        libMesh::Real k = this->_k(T);
        libMesh::Real cp = this->_cp(T);

        libMesh::Number rho_cp = rho*this->_cp(T);

        libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u());

        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::Real tau_E = this->_stab_helper.compute_tau_energy( context, qp, g, G, rho, U, k, cp, this->_is_steady );

        libMesh::Real RE_s = this->compute_res_energy_steady( context, qp );

        for (unsigned int i=0; i != n_T_dofs; i++)
          {
            FT(i) += ( rho_cp*tau_E*RE_s*U*T_gradphi[i][qp]
                       + tau_E*RE_s*k*(T_hessphi[i][qp](0,0) + T_hessphi[i][qp](1,1) + T_hessphi[i][qp](2,2) )
                       )*JxW[qp];
          }

      }

    return;
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokesBraackStabilization<Mu,SH,TC>::assemble_continuity_mass_residual( bool /*compute_jacobian*/,
                                                                                            AssemblyContext& context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.get_dof_indices(this->_press_var.p()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
      context.get_element_fe(this->_press_var.p())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(this->_press_var.p()); // R_{p}

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u());

        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::Real T = context.fixed_interior_value( this->_temp_vars.T(), qp );
        libMesh::Real rho = this->rho( T, this->get_p0_transient( context, qp ) );

        libMesh::Real mu = this->_mu(T);
        libMesh::Real k = this->_k(T);
        libMesh::Real cp = this->_cp(T);

        libMesh::RealGradient U( context.fixed_interior_value( this->_flow_vars.u(), qp ),
                                 context.fixed_interior_value( this->_flow_vars.v(), qp ) );
        if( this->_flow_vars.dim() == 3 )
          U(2) = context.fixed_interior_value( this->_flow_vars.w(), qp );

        libMesh::Real tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, rho, U, mu, false );
        libMesh::RealGradient RM_t = this->compute_res_momentum_transient( context, qp );

        libMesh::Real tau_E = this->_stab_helper.compute_tau_energy( context, qp, g, G, rho, U, k, cp, false );
        libMesh::Real RE_t = this->compute_res_energy_transient( context, qp );

        // Now a loop over the pressure degrees of freedom.  This
        // computes the contributions of the continuity equation.
        for (unsigned int i=0; i != n_p_dofs; i++)
          {
            Fp(i) += ( tau_M*RM_t*p_dphi[i][qp]
                       +  tau_E*RE_t*(U*p_dphi[i][qp])/T
                       )*JxW[qp];
          }
      }

    return;
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokesBraackStabilization<Mu,SH,TC>::assemble_momentum_mass_residual( bool /*compute_jacobian*/,
                                                                                          AssemblyContext& context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();

    // Check number of dofs is same for _flow_vars.u(), v_var and w_var.
    libmesh_assert (n_u_dofs == context.get_dof_indices(this->_flow_vars.v()).size());
    if (this->_flow_vars.dim() == 3)
      libmesh_assert (n_u_dofs == context.get_dof_indices(this->_flow_vars.w()).size());

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // The velocity shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_flow_vars.u())->get_dphi();

    const std::vector<std::vector<libMesh::RealTensor> >& u_hessphi =
      context.get_element_fe(this->_flow_vars.u())->get_d2phi();

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u()); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v()); // R_{v}
    libMesh::DenseSubVector<libMesh::Real>* Fw = NULL;

    if( this->_flow_vars.dim() == 3 )
      {
        Fw  = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Real T = context.fixed_interior_value( this->_temp_vars.T(), qp );
        libMesh::Real rho = this->rho( T, this->get_p0_transient( context, qp ) );

        libMesh::Real mu = this->_mu(T);

        libMesh::RealGradient U( context.fixed_interior_value(this->_flow_vars.u(), qp),
                                 context.fixed_interior_value(this->_flow_vars.v(), qp) );

        libMesh::RealGradient grad_u = context.fixed_interior_gradient(this->_flow_vars.u(), qp);
        libMesh::RealGradient grad_v = context.fixed_interior_gradient(this->_flow_vars.v(), qp);
        libMesh::RealGradient grad_w;

        if( this->_flow_vars.dim() == 3 )
          {
            U(2) = context.fixed_interior_value(this->_flow_vars.w(), qp);
            grad_w = context.fixed_interior_gradient(this->_flow_vars.w(), qp);
          }

        libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u());

        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::Real tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, rho, U, mu, false );
        libMesh::Real tau_C = this->_stab_helper.compute_tau_continuity( tau_M, g );

        libMesh::Real RC_t = this->compute_res_continuity_transient( context, qp );
        libMesh::RealGradient RM_t = this->compute_res_momentum_transient( context, qp );

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) += (  tau_C*RC_t*u_gradphi[i][qp](0)
                        + tau_M*RM_t(0)*rho*U*u_gradphi[i][qp]
                        + mu*tau_M*RM_t(0)*(u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1)
                                            + u_hessphi[i][qp](0,0) + u_hessphi[i][qp](0,1)
                                            - 2.0/3.0*(u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,0)) )
                        )*JxW[qp];

            Fv(i) += ( tau_C*RC_t*u_gradphi[i][qp](1)
                       + tau_M*RM_t(1)*rho*U*u_gradphi[i][qp]
                       + mu*tau_M*RM_t(1)*(u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1)
                                           + u_hessphi[i][qp](1,0) + u_hessphi[i][qp](1,1)
                                           - 2.0/3.0*(u_hessphi[i][qp](0,1) + u_hessphi[i][qp](1,1)) )
                       )*JxW[qp];

            if( this->_flow_vars.dim() == 3 )
              {
                Fu(i) -= mu*tau_M*RM_t(0)*(u_hessphi[i][qp](2,2) + u_hessphi[i][qp](0,2)
                                           - 2.0/3.0*u_hessphi[i][qp](2,0))*JxW[qp];

                Fv(i) -= mu*tau_M*RM_t(1)*(u_hessphi[i][qp](2,2) + u_hessphi[i][qp](1,2)
                                           - 2.0/3.0*u_hessphi[i][qp](2,1))*JxW[qp];

                (*Fw)(i) -= ( tau_C*RC_t*u_gradphi[i][qp](2)
                              + tau_M*RM_t(2)*rho*U*u_gradphi[i][qp]
                              + mu*tau_M*RM_t(2)*(u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) + u_hessphi[i][qp](2,2)
                                                  + u_hessphi[i][qp](2,0) + u_hessphi[i][qp](2,1) + u_hessphi[i][qp](2,2)
                                                  - 2.0/3.0*(u_hessphi[i][qp](0,2) + u_hessphi[i][qp](1,2)
                                                             + u_hessphi[i][qp](2,2))
                                                  )
                              )*JxW[qp];
              }
          }

      }
    return;
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokesBraackStabilization<Mu,SH,TC>::assemble_energy_mass_residual( bool /*compute_jacobian*/,
                                                                                        AssemblyContext& context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_temp_vars.T())->get_JxW();

    // The temperature shape functions gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(this->_temp_vars.T())->get_dphi();

    const std::vector<std::vector<libMesh::RealTensor> >& T_hessphi =
      context.get_element_fe(this->_temp_vars.T())->get_d2phi();

    libMesh::DenseSubVector<libMesh::Number> &FT = context.get_elem_residual(this->_temp_vars.T()); // R_{T}

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Number u, v;
        u = context.fixed_interior_value(this->_flow_vars.u(), qp);
        v = context.fixed_interior_value(this->_flow_vars.v(), qp);

        libMesh::Gradient grad_T = context.fixed_interior_gradient(this->_temp_vars.T(), qp);

        libMesh::NumberVectorValue U(u,v);
        if (this->_flow_vars.dim() == 3)
          U(2) = context.fixed_interior_value(this->_flow_vars.w(), qp); // w

        libMesh::Real T = context.fixed_interior_value( this->_temp_vars.T(), qp );
        libMesh::Real rho = this->rho( T, this->get_p0_transient( context, qp ) );

        libMesh::Real k = this->_k(T);
        libMesh::Real cp = this->_cp(T);

        libMesh::Number rho_cp = rho*cp;

        libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u());

        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::Real tau_E = this->_stab_helper.compute_tau_energy( context, qp, g, G, rho, U, k, cp, false );

        libMesh::Real RE_t = this->compute_res_energy_transient( context, qp );

        for (unsigned int i=0; i != n_T_dofs; i++)
          {
            FT(i) += ( rho_cp*tau_E*RE_t*U*T_gradphi[i][qp]
                       + tau_E*RE_t*k*(T_hessphi[i][qp](0,0) + T_hessphi[i][qp](1,1) + T_hessphi[i][qp](2,2))
                       )*JxW[qp];
          }

      }

    return;
  }


}  // namespace GRINS

// Instantiate
template class GRINS::LowMachNavierStokesBraackStabilization<GRINS::ConstantViscosity,GRINS::ConstantSpecificHeat,GRINS::ConstantConductivity>;
