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
#include "grins/boussinesq_buoyancy_adjoint_stab.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/inc_nav_stokes_macro.h"
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace GRINS
{
  template<class Mu>
  BoussinesqBuoyancyAdjointStabilization<Mu>::BoussinesqBuoyancyAdjointStabilization( const std::string& physics_name, const GetPot& input )
    : BoussinesqBuoyancyBase(physics_name,input),
      _mu(input,MaterialsParsing::material_name(input,PhysicsNaming::boussinesq_buoyancy())),
      _stab_helper( physics_name+"StabHelper", input )
  {}

  template<class Mu>
  BoussinesqBuoyancyAdjointStabilization<Mu>::~BoussinesqBuoyancyAdjointStabilization()
  {
    return;
  }

  template<class Mu>
  void BoussinesqBuoyancyAdjointStabilization<Mu>::init_context( AssemblyContext & context )
  {
    context.get_element_fe(this->_press_var.p())->get_dphi();

    context.get_element_fe(this->_flow_vars.u())->get_dphi();
    context.get_element_fe(this->_flow_vars.u())->get_d2phi();

    return;
  }

  template<class Mu>
  void BoussinesqBuoyancyAdjointStabilization<Mu>::element_time_derivative
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.get_dof_indices(_flow_vars.u()).size();
    const unsigned int n_T_dofs = context.get_dof_indices(_temp_vars.T()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(_flow_vars.u())->get_JxW();

    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_temp_vars.T())->get_phi();

    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_flow_vars.u())->get_dphi();

    const std::vector<std::vector<libMesh::RealTensor> >& u_hessphi =
      context.get_element_fe(this->_flow_vars.u())->get_d2phi();

    // Get residuals and jacobians
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(_flow_vars.u()); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(_flow_vars.v()); // R_{v}
    libMesh::DenseSubVector<libMesh::Number> *Fw = NULL;

    libMesh::DenseSubMatrix<libMesh::Number> &KuT =
      context.get_elem_jacobian(_flow_vars.u(), _temp_vars.T()); // J_{uT}
    libMesh::DenseSubMatrix<libMesh::Number> &KvT =
      context.get_elem_jacobian(_flow_vars.v(), _temp_vars.T()); // J_{vT}
    libMesh::DenseSubMatrix<libMesh::Number> &Kuu =
      context.get_elem_jacobian(_flow_vars.u(), _flow_vars.u()); // J_{uu}
    libMesh::DenseSubMatrix<libMesh::Number> &Kuv =
      context.get_elem_jacobian(_flow_vars.u(), _flow_vars.v()); // J_{uv}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvu =
      context.get_elem_jacobian(_flow_vars.v(), _flow_vars.u()); // J_{vu}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvv =
      context.get_elem_jacobian(_flow_vars.v(), _flow_vars.v()); // J_{vv}

    libMesh::DenseSubMatrix<libMesh::Number> *KwT = NULL;
    libMesh::DenseSubMatrix<libMesh::Number> *Kuw = NULL;
    libMesh::DenseSubMatrix<libMesh::Number> *Kvw = NULL;
    libMesh::DenseSubMatrix<libMesh::Number> *Kwu = NULL;
    libMesh::DenseSubMatrix<libMesh::Number> *Kwv = NULL;
    libMesh::DenseSubMatrix<libMesh::Number> *Kww = NULL;


    if(this->_flow_vars.dim() == 3)
      {
        Fw = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}
        KwT = &context.get_elem_jacobian
          (_flow_vars.w(), _temp_vars.T()); // J_{wT}
        Kuw = &context.get_elem_jacobian
          (_flow_vars.u(), _flow_vars.w()); // J_{uw}
        Kvw = &context.get_elem_jacobian
          (_flow_vars.v(), _flow_vars.w()); // J_{vw}
        Kwu = &context.get_elem_jacobian
          (_flow_vars.w(), _flow_vars.u()); // J_{wu}
        Kwv = &context.get_elem_jacobian
          (_flow_vars.w(), _flow_vars.v()); // J_{wv}
        Kww = &context.get_elem_jacobian
          (_flow_vars.w(), _flow_vars.w()); // J_{ww}
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

        // Compute the solution & its gradient at the old Newton iterate.
        libMesh::Number T;
        T = context.interior_value(_temp_vars.T(), qp);

        libMesh::RealGradient d_residual_dT = _rho*_beta_T*_g;
        // d_residual_dU = 0
        libMesh::RealGradient residual = (T-_T_ref)*d_residual_dT;

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            libMesh::Real test_func = this->_rho*U*u_gradphi[i][qp] +
              mu_qp*( u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) + u_hessphi[i][qp](2,2) );
            Fu(i) += -tau_M*residual(0)*test_func*JxW[qp];

            Fv(i) += -tau_M*residual(1)*test_func*JxW[qp];

            if (this->_flow_vars.dim() == 3)
              {
                (*Fw)(i) += -tau_M*residual(2)*test_func*JxW[qp];
              }

            if (compute_jacobian)
              {
                libMesh::Gradient d_test_func_dU = this->_rho*u_gradphi[i][qp];
                // d_test_func_dT = 0

                for (unsigned int j=0; j != n_u_dofs; ++j)
                  {
                    Kuu(i,j) += -tau_M*residual(0)*d_test_func_dU(0)*u_phi[j][qp]*JxW[qp] * context.get_elem_solution_derivative();
                    Kuu(i,j) += -d_tau_M_dU(0)*u_phi[j][qp]*residual(0)*test_func*JxW[qp] * context.get_elem_solution_derivative();
                    Kuv(i,j) += -tau_M*residual(0)*d_test_func_dU(1)*u_phi[j][qp]*JxW[qp] * context.get_elem_solution_derivative();
                    Kuv(i,j) += -d_tau_M_dU(1)*u_phi[j][qp]*residual(0)*test_func*JxW[qp] * context.get_elem_solution_derivative();
                    Kvu(i,j) += -tau_M*residual(1)*d_test_func_dU(0)*u_phi[j][qp]*JxW[qp] * context.get_elem_solution_derivative();
                    Kvu(i,j) += -d_tau_M_dU(0)*u_phi[j][qp]*residual(1)*test_func*JxW[qp] * context.get_elem_solution_derivative();
                    Kvv(i,j) += -tau_M*residual(1)*d_test_func_dU(1)*u_phi[j][qp]*JxW[qp] * context.get_elem_solution_derivative();
                    Kvv(i,j) += -d_tau_M_dU(1)*u_phi[j][qp]*residual(1)*test_func*JxW[qp] * context.get_elem_solution_derivative();
                  }

                for (unsigned int j=0; j != n_T_dofs; ++j)
                  {
                    // KuT(i,j) += -tau_M*residual(0)*dtest_func_dT[j]*JxW[qp] * context.get_elem_solution_derivative();
                    KuT(i,j) += -tau_M*d_residual_dT(0)*T_phi[j][qp]*test_func*JxW[qp] * context.get_elem_solution_derivative();
                    // KvT(i,j) += -tau_M*residual(1)*dtest_func_dT[j]*JxW[qp] * context.get_elem_solution_derivative();
                    KvT(i,j) += -tau_M*d_residual_dT(1)*T_phi[j][qp]*test_func*JxW[qp] * context.get_elem_solution_derivative();
                  }
                if (this->_flow_vars.dim() == 3)
                  {
                    for (unsigned int j=0; j != n_T_dofs; ++j)
                      {
                        // KwT(i,j) += -tau_M*residual(2)*dtest_func_dT[j]*JxW[qp] * context.get_elem_solution_derivative();
                        (*KwT)(i,j) += -tau_M*d_residual_dT(2)*T_phi[j][qp]*test_func*JxW[qp] * context.get_elem_solution_derivative();
                      }

                    for (unsigned int j=0; j != n_u_dofs; ++j)
                      {
                        (*Kuw)(i,j) += -tau_M*residual(0)*d_test_func_dU(2)*u_phi[j][qp]*JxW[qp] * context.get_elem_solution_derivative();
                        (*Kuw)(i,j) += -d_tau_M_dU(2)*u_phi[j][qp]*residual(0)*test_func*JxW[qp] * context.get_elem_solution_derivative();
                        (*Kvw)(i,j) += -tau_M*residual(1)*d_test_func_dU(2)*u_phi[j][qp]*JxW[qp] * context.get_elem_solution_derivative();
                        (*Kvw)(i,j) += -d_tau_M_dU(2)*u_phi[j][qp]*residual(1)*test_func*JxW[qp] * context.get_elem_solution_derivative();
                        (*Kwu)(i,j) += -tau_M*residual(2)*d_test_func_dU(0)*u_phi[j][qp]*JxW[qp] * context.get_elem_solution_derivative();
                        (*Kwu)(i,j) += -d_tau_M_dU(0)*u_phi[j][qp]*residual(2)*test_func*JxW[qp] * context.get_elem_solution_derivative();
                        (*Kwv)(i,j) += -tau_M*residual(2)*d_test_func_dU(1)*u_phi[j][qp]*JxW[qp] * context.get_elem_solution_derivative();
                        (*Kwv)(i,j) += -d_tau_M_dU(1)*u_phi[j][qp]*residual(2)*test_func*JxW[qp] * context.get_elem_solution_derivative();
                        (*Kww)(i,j) += -tau_M*residual(2)*d_test_func_dU(2)*u_phi[j][qp]*JxW[qp] * context.get_elem_solution_derivative();
                        (*Kww)(i,j) += -d_tau_M_dU(2)*u_phi[j][qp]*residual(2)*test_func*JxW[qp] * context.get_elem_solution_derivative();
                      }
                  }

              } // End compute_jacobian check

          } // End i dof loop
      } // End quadrature loop
  }

  template<class Mu>
  void BoussinesqBuoyancyAdjointStabilization<Mu>::element_constraint
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.get_dof_indices(_press_var.p()).size();
    const unsigned int n_u_dofs = context.get_dof_indices(_flow_vars.u()).size();
    const unsigned int n_T_dofs = context.get_dof_indices(_temp_vars.T()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(_flow_vars.u())->get_JxW();

    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_temp_vars.T())->get_phi();

    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
      context.get_element_fe(this->_press_var.p())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(this->_press_var.p()); // R_{p}

    libMesh::DenseSubMatrix<libMesh::Number> &KpT =
      context.get_elem_jacobian(_press_var.p(), _temp_vars.T()); // J_{pT}
    libMesh::DenseSubMatrix<libMesh::Number> &Kpu =
      context.get_elem_jacobian(_press_var.p(), _flow_vars.u()); // J_{pu}
    libMesh::DenseSubMatrix<libMesh::Number> &Kpv =
      context.get_elem_jacobian(_press_var.p(), _flow_vars.v()); // J_{pv}
    libMesh::DenseSubMatrix<libMesh::Number> *Kpw = NULL;

    if(this->_flow_vars.dim() == 3)
      {
        Kpw = &context.get_elem_jacobian
          (_press_var.p(), _flow_vars.w()); // J_{pw}
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

        // Compute the solution & its gradient at the old Newton iterate.
        libMesh::Number T;
        T = context.interior_value(_temp_vars.T(), qp);

        libMesh::RealGradient d_residual_dT = _rho*_beta_T*_g;
        // d_residual_dU = 0
        libMesh::RealGradient residual = (T-_T_ref)*d_residual_dT;

        // First, an i-loop over the velocity degrees of freedom.
        // We know that n_u_dofs == n_v_dofs so we can compute contributions
        // for both at the same time.
        for (unsigned int i=0; i != n_p_dofs; i++)
          {
            Fp(i) += tau_M*residual*p_dphi[i][qp]*JxW[qp];

            if (compute_jacobian)
              {
                for (unsigned int j=0; j != n_T_dofs; ++j)
                  {
                    KpT(i,j) += tau_M*d_residual_dT*T_phi[j][qp]*p_dphi[i][qp]*JxW[qp] * context.get_elem_solution_derivative();
                  }

                for (unsigned int j=0; j != n_u_dofs; ++j)
                  {
                    Kpu(i,j) += d_tau_M_dU(0)*u_phi[j][qp]*residual*p_dphi[i][qp]*JxW[qp] * context.get_elem_solution_derivative();
                    Kpv(i,j) += d_tau_M_dU(1)*u_phi[j][qp]*residual*p_dphi[i][qp]*JxW[qp] * context.get_elem_solution_derivative();
                  }
                if( this->_flow_vars.dim() == 3 )
                  for (unsigned int j=0; j != n_u_dofs; ++j)
                    {
                      (*Kpw)(i,j) += d_tau_M_dU(2)*u_phi[j][qp]*residual*p_dphi[i][qp]*JxW[qp] * context.get_elem_solution_derivative();
                    }
              }
          }
      } // End quadrature loop
  }

  template<class Mu>
  void BoussinesqBuoyancyAdjointStabilization<Mu>::register_parameter
  ( const std::string & param_name,
    libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
    const
  {
    ParameterUser::register_parameter(param_name, param_pointer);
    _mu.register_parameter(param_name, param_pointer);
  }


} // namespace GRINS

// Instantiate
INSTANTIATE_INC_NS_SUBCLASS(BoussinesqBuoyancyAdjointStabilization);
