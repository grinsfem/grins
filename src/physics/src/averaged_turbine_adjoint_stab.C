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
#include "grins_config.h"
#include "grins/averaged_turbine_adjoint_stab.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/inc_nav_stokes_macro.h"

// libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{

  template<class Mu>
  AveragedTurbineAdjointStabilization<Mu>::AveragedTurbineAdjointStabilization( const std::string& physics_name, const GetPot& input )
    : AveragedTurbineBase<Mu>(physics_name, input),
    _stab_helper( physics_name+"StabHelper", input )
  {}

  template<class Mu>
  AveragedTurbineAdjointStabilization<Mu>::~AveragedTurbineAdjointStabilization()
  {
    return;
  }

  template<class Mu>
  void AveragedTurbineAdjointStabilization<Mu>::init_context( AssemblyContext& context )
  {
    context.get_element_fe(this->_press_var.p())->get_dphi();

    context.get_element_fe(this->_flow_vars.u())->get_xyz();
    context.get_element_fe(this->_flow_vars.u())->get_phi();
    context.get_element_fe(this->_flow_vars.u())->get_dphi();
    context.get_element_fe(this->_flow_vars.u())->get_d2phi();

    return;
  }

  template<class Mu>
  void AveragedTurbineAdjointStabilization<Mu>::element_time_derivative
  ( bool compute_jacobian,
    AssemblyContext & context )
  {
    libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u());

    // Element Jacobian * quadrature weights for interior integration
    const std::vector<libMesh::Real> &JxW = fe->get_JxW();

    // The shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& u_phi = fe->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_flow_vars.u())->get_dphi();

    const std::vector<std::vector<libMesh::RealTensor> >& u_hessphi =
      context.get_element_fe(this->_flow_vars.u())->get_d2phi();

    const std::vector<libMesh::Point>& u_qpoint = fe->get_xyz();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubMatrix<libMesh::Number> &Kuu = context.get_elem_jacobian(this->_flow_vars.u(), this->_flow_vars.u()); // R_{u},{u}
    libMesh::DenseSubMatrix<libMesh::Number> &Kuv = context.get_elem_jacobian(this->_flow_vars.u(), this->_flow_vars.v()); // R_{u},{v}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvu = context.get_elem_jacobian(this->_flow_vars.v(), this->_flow_vars.u()); // R_{v},{u}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvv = context.get_elem_jacobian(this->_flow_vars.v(), this->_flow_vars.v()); // R_{v},{v}

    libMesh::DenseSubMatrix<libMesh::Number> &Kus =
      context.get_elem_jacobian(this->_flow_vars.u(),
                                this->fan_speed_var()); // R_{u},{s}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvs =
      context.get_elem_jacobian(this->_flow_vars.v(),
                                this->fan_speed_var()); // R_{v},{s}

    libMesh::DenseSubMatrix<libMesh::Number>* Kwu = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kwv = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kww = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kuw = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kvw = NULL;

    libMesh::DenseSubMatrix<libMesh::Number>* Kws = NULL;

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u()); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v()); // R_{v}
    libMesh::DenseSubVector<libMesh::Number>* Fw = NULL;

    if( this->_flow_vars.dim() == 3 )
      {
        Kuw = &context.get_elem_jacobian(this->_flow_vars.u(), this->_flow_vars.w()); // R_{u},{w}
        Kvw = &context.get_elem_jacobian(this->_flow_vars.v(), this->_flow_vars.w()); // R_{v},{w}

        Kwu = &context.get_elem_jacobian(this->_flow_vars.w(), this->_flow_vars.u()); // R_{w},{u}
        Kwv = &context.get_elem_jacobian(this->_flow_vars.w(), this->_flow_vars.v()); // R_{w},{v}
        Kww = &context.get_elem_jacobian(this->_flow_vars.w(), this->_flow_vars.w()); // R_{w},{w}

        Kws = &context.get_elem_jacobian(this->_flow_vars.w(), this->fan_speed_var()); // R_{w},{s}

        Fw  = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();

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

        libMesh::Number s = context.interior_value(this->fan_speed_var(), qp);

        // Compute the viscosity at this qp
        libMesh::Real mu_qp = this->_mu(context, qp);

        libMesh::Real tau_M;
        libMesh::Real d_tau_M_d_rho;
        libMesh::Gradient d_tau_M_dU;

        if (compute_jacobian)
          this->_stab_helper.compute_tau_momentum_and_derivs
            ( context, qp, g, G, this->_rho, U, mu_qp,
              tau_M, d_tau_M_d_rho, d_tau_M_dU,
              this->_is_steady );
        else
          tau_M = this->_stab_helper.compute_tau_momentum
            ( context, qp, g, G, this->_rho, U, mu_qp,
              this->_is_steady );

        libMesh::NumberVectorValue U_B_1;
        libMesh::NumberVectorValue F;
        libMesh::NumberTensorValue dFdU;
        libMesh::NumberTensorValue* dFdU_ptr =
          compute_jacobian ? &dFdU : NULL;
        libMesh::NumberVectorValue dFds;
        libMesh::NumberVectorValue* dFds_ptr =
          compute_jacobian ? &dFds : NULL;
        if (!this->compute_force(u_qpoint[qp], context.time, U, s,
                                 U_B_1, F, dFdU_ptr, dFds_ptr))
          continue;

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            libMesh::Real test_func = this->_rho*U*u_gradphi[i][qp] +
              mu_qp*( u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) + u_hessphi[i][qp](2,2) );
            Fu(i) += tau_M*F(0)*test_func*JxW[qp];

            Fv(i) += tau_M*F(1)*test_func*JxW[qp];

            if (this->_flow_vars.dim() == 3)
              {
                (*Fw)(i) += tau_M*F(2)*test_func*JxW[qp];
              }

            if (compute_jacobian)
              {
                libMesh::Gradient d_test_func_dU = this->_rho*u_gradphi[i][qp];

                for (unsigned int j=0; j != n_u_dofs; ++j)
                  {
                    Kuu(i,j) += tau_M*F(0)*d_test_func_dU(0)*u_phi[j][qp]*JxW[qp];
                    Kuu(i,j) += d_tau_M_dU(0)*u_phi[j][qp]*F(0)*test_func*JxW[qp];
                    Kuu(i,j) += tau_M*dFdU(0,0)*u_phi[j][qp]*test_func*JxW[qp];
                    Kuv(i,j) += tau_M*F(0)*d_test_func_dU(1)*u_phi[j][qp]*JxW[qp];
                    Kuv(i,j) += d_tau_M_dU(1)*u_phi[j][qp]*F(0)*test_func*JxW[qp];
                    Kuv(i,j) += tau_M*dFdU(0,1)*u_phi[j][qp]*test_func*JxW[qp];
                    Kvu(i,j) += tau_M*F(1)*d_test_func_dU(0)*u_phi[j][qp]*JxW[qp];
                    Kvu(i,j) += d_tau_M_dU(0)*u_phi[j][qp]*F(1)*test_func*JxW[qp];
                    Kvu(i,j) += tau_M*dFdU(1,0)*u_phi[j][qp]*test_func*JxW[qp];
                    Kvv(i,j) += tau_M*F(1)*d_test_func_dU(1)*u_phi[j][qp]*JxW[qp];
                    Kvv(i,j) += d_tau_M_dU(1)*u_phi[j][qp]*F(1)*test_func*JxW[qp];
                    Kvv(i,j) += tau_M*dFdU(1,1)*u_phi[j][qp]*test_func*JxW[qp];

                    Kus(i,0)   += tau_M*dFds(0)*test_func*JxW[qp];
                    Kvs(i,0)   += tau_M*dFds(1)*test_func*JxW[qp];
                  }

                if (this->_flow_vars.dim() == 3)
                  {
                    for (unsigned int j=0; j != n_u_dofs; ++j)
                      {
                        (*Kuw)(i,j) += tau_M*F(0)*d_test_func_dU(2)*u_phi[j][qp]*JxW[qp];
                        (*Kuw)(i,j) += d_tau_M_dU(2)*u_phi[j][qp]*F(0)*test_func*JxW[qp];
                        (*Kuw)(i,j) += tau_M*dFdU(0,2)*u_phi[j][qp]*test_func*JxW[qp];
                        (*Kvw)(i,j) += tau_M*F(1)*d_test_func_dU(2)*u_phi[j][qp]*JxW[qp];
                        (*Kvw)(i,j) += d_tau_M_dU(2)*u_phi[j][qp]*F(1)*test_func*JxW[qp];
                        (*Kvw)(i,j) += tau_M*dFdU(1,2)*u_phi[j][qp]*test_func*JxW[qp];
                        (*Kwu)(i,j) += tau_M*F(2)*d_test_func_dU(0)*u_phi[j][qp]*JxW[qp];
                        (*Kwu)(i,j) += d_tau_M_dU(0)*u_phi[j][qp]*F(2)*test_func*JxW[qp];
                        (*Kwu)(i,j) += tau_M*dFdU(2,0)*u_phi[j][qp]*test_func*JxW[qp];
                        (*Kwv)(i,j) += tau_M*F(2)*d_test_func_dU(1)*u_phi[j][qp]*JxW[qp];
                        (*Kwv)(i,j) += d_tau_M_dU(1)*u_phi[j][qp]*F(2)*test_func*JxW[qp];
                        (*Kwv)(i,j) += tau_M*dFdU(2,1)*u_phi[j][qp]*test_func*JxW[qp];
                        (*Kww)(i,j) += tau_M*F(2)*d_test_func_dU(2)*u_phi[j][qp]*JxW[qp];
                        (*Kww)(i,j) += d_tau_M_dU(2)*u_phi[j][qp]*F(2)*test_func*JxW[qp];
                        (*Kww)(i,j) += tau_M*dFdU(2,2)*u_phi[j][qp]*test_func*JxW[qp];

                        (*Kws)(i,0) += tau_M*dFds(2)*test_func*JxW[qp];
                      }
                  }

              } // End compute_jacobian check

          } // End i dof loop
      }
  }

  template<class Mu>
  void AveragedTurbineAdjointStabilization<Mu>::element_constraint
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.get_dof_indices(this->_press_var.p()).size();
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    const std::vector<libMesh::Point>& u_qpoint =
      context.get_element_fe(this->_flow_vars.u())->get_xyz();

    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
      context.get_element_fe(this->_press_var.p())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(this->_press_var.p()); // R_{p}

    libMesh::DenseSubMatrix<libMesh::Number> &Kpu =
      context.get_elem_jacobian(this->_press_var.p(), this->_flow_vars.u()); // J_{pu}
    libMesh::DenseSubMatrix<libMesh::Number> &Kpv =
      context.get_elem_jacobian(this->_press_var.p(), this->_flow_vars.v()); // J_{pv}
    libMesh::DenseSubMatrix<libMesh::Number> &Kps =
      context.get_elem_jacobian(this->_press_var.p(), this->fan_speed_var()); // J_{ps}

    libMesh::DenseSubMatrix<libMesh::Number> *Kpw = NULL;

    if(this->_flow_vars.dim() == 3)
      {
        Kpw = &context.get_elem_jacobian
          (this->_press_var.p(), this->_flow_vars.w()); // J_{pw}
      }

    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
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

        libMesh::Number s = context.interior_value(this->fan_speed_var(), qp);

        // Compute the viscosity at this qp
        libMesh::Real mu_qp = this->_mu(context, qp);

        libMesh::Real tau_M;
        libMesh::Real d_tau_M_d_rho;
        libMesh::Gradient d_tau_M_dU;

        if (compute_jacobian)
          this->_stab_helper.compute_tau_momentum_and_derivs
            ( context, qp, g, G, this->_rho, U, mu_qp,
              tau_M, d_tau_M_d_rho, d_tau_M_dU,
              this->_is_steady );
        else
          tau_M = this->_stab_helper.compute_tau_momentum
            ( context, qp, g, G, this->_rho, U, mu_qp,
              this->_is_steady );

        libMesh::NumberVectorValue U_B_1;
        libMesh::NumberVectorValue F;
        libMesh::NumberTensorValue dFdU;
        libMesh::NumberTensorValue* dFdU_ptr =
          compute_jacobian ? &dFdU : NULL;
        libMesh::NumberVectorValue dFds;
        libMesh::NumberVectorValue* dFds_ptr =
          compute_jacobian ? &dFds : NULL;
        if (!this->compute_force(u_qpoint[qp], context.time, U, s,
                                 U_B_1, F, dFdU_ptr, dFds_ptr))
          continue;

        // First, an i-loop over the velocity degrees of freedom.
        // We know that n_u_dofs == n_v_dofs so we can compute contributions
        // for both at the same time.
        for (unsigned int i=0; i != n_p_dofs; i++)
          {
            Fp(i) += -tau_M*(F*p_dphi[i][qp])*JxW[qp];

            if (compute_jacobian)
              {
                Kps(i,0) += -tau_M*(dFds*p_dphi[i][qp])*JxW[qp];

                for (unsigned int j=0; j != n_u_dofs; ++j)
                  {
                    Kpu(i,j) += -d_tau_M_dU(0)*u_phi[j][qp]*(F*p_dphi[i][qp])*JxW[qp];
                    Kpv(i,j) += -d_tau_M_dU(1)*u_phi[j][qp]*(F*p_dphi[i][qp])*JxW[qp];
                    for (unsigned int d=0; d != 3; ++d)
                      {
                        Kpu(i,j) += -tau_M*dFdU(d,0)*u_phi[j][qp]*p_dphi[i][qp](d)*JxW[qp];
                        Kpv(i,j) += -tau_M*dFdU(d,1)*u_phi[j][qp]*p_dphi[i][qp](d)*JxW[qp];
                      }
                  }
                if( this->_flow_vars.dim() == 3 )
                  for (unsigned int j=0; j != n_u_dofs; ++j)
                    {
                      (*Kpw)(i,j) += -d_tau_M_dU(2)*u_phi[j][qp]*F*p_dphi[i][qp]*JxW[qp];
                      for (unsigned int d=0; d != 3; ++d)
                        {
                          (*Kpw)(i,j) += -tau_M*dFdU(d,2)*u_phi[j][qp]*p_dphi[i][qp](d)*JxW[qp];
                        }
                    }
              }
          }
      } // End quadrature loop
  }


} // namespace GRINS

// Instantiate
INSTANTIATE_INC_NS_SUBCLASS(AveragedTurbineAdjointStabilization);
