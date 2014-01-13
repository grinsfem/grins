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
#include "grins/assembly_context.h"
#include "grins/heat_transfer_adjoint_stab.h"

// libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{

  HeatTransferAdjointStabilization::HeatTransferAdjointStabilization( const std::string& physics_name, 
                                                                      const GetPot& input )
    : HeatTransferStabilizationBase(physics_name,input)
  {
    return;
  }

  HeatTransferAdjointStabilization::~HeatTransferAdjointStabilization()
  {
    return;
  }

  void HeatTransferAdjointStabilization::element_time_derivative( bool compute_jacobian,
                                                                  AssemblyContext& context,
                                                                  CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("HeatTransferAdjointStabilization::element_time_derivative");
#endif

    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T_var()).size();
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u_var()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_temp_vars.T_var())->get_JxW();

    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u_var())->get_phi();

    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_temp_vars.T_var())->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(this->_temp_vars.T_var())->get_dphi();

    const std::vector<std::vector<libMesh::RealTensor> >& T_hessphi =
      context.get_element_fe(this->_temp_vars.T_var())->get_d2phi();

    /*
      const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u_var()).size();

      const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u_var())->get_phi();

      libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u_var()); // R_{p}
      libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v_var()); // R_{p}
      libMesh::DenseSubVector<libMesh::Number> *Fw = NULL;
      if(this->_dim == 3)
      {
      Fw = &context.get_elem_residual(this->_flow_vars.w_var()); // R_{w}
      }
    */

    libMesh::DenseSubVector<libMesh::Number> &FT = context.get_elem_residual(this->_temp_vars.T_var()); // R_{T}
    libMesh::DenseSubMatrix<libMesh::Number> &KTT = 
      context.get_elem_jacobian(_temp_vars.T_var(), _temp_vars.T_var()); // J_{TT}
    libMesh::DenseSubMatrix<libMesh::Number> &KTu = 
      context.get_elem_jacobian(_temp_vars.T_var(), _flow_vars.u_var()); // J_{Tu}
    libMesh::DenseSubMatrix<libMesh::Number> &KTv = 
      context.get_elem_jacobian(_temp_vars.T_var(), _flow_vars.v_var()); // J_{Tv}
    libMesh::DenseSubMatrix<libMesh::Number> *KTw = NULL;

    if(this->_dim == 3)
      {
        KTw = &context.get_elem_jacobian
          (_temp_vars.T_var(), _flow_vars.w_var()); // J_{Tw}
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::FEBase* fe = context.get_element_fe(this->_temp_vars.T_var());

        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::RealGradient U( context.interior_value( this->_flow_vars.u_var(), qp ),
                                 context.interior_value( this->_flow_vars.v_var(), qp ) );
        if( this->_dim == 3 )
          {
            U(2) = context.interior_value( this->_flow_vars.w_var(), qp );
          }
      
        //libMesh::RealGradient grad_T = context.interior_gradient( this->_temp_vars.T_var(), qp );

        libMesh::Real tau_E, RE_s;
        libMesh::Real d_tau_E_d_rho, d_RE_s_dT;
        libMesh::Gradient d_tau_E_dU, d_RE_s_dgradT, d_RE_s_dU;
        libMesh::Tensor d_RE_s_dhessT;
       
        if (compute_jacobian)
          {
            this->_stab_helper.compute_tau_energy_and_derivs
              ( context, G, _rho, _Cp, _k,  U,
                tau_E, d_tau_E_d_rho, d_tau_E_dU, this->_is_steady );
            this->_stab_helper.compute_res_energy_steady_and_derivs
              ( context, qp, _rho, _Cp, _k,
                RE_s, d_RE_s_dT, d_RE_s_dgradT, d_RE_s_dhessT,
                d_RE_s_dU );
          }
        else
          {
            tau_E = this->_stab_helper.compute_tau_energy
              ( context, G, _rho, _Cp, _k,  U, this->_is_steady );

            RE_s = this->_stab_helper.compute_res_energy_steady
              ( context, qp, _rho, _Cp, _k );
          }

        /*
          for (unsigned int i=0; i != n_u_dofs; i++)
          {
          Fu(i) += -tau_E*RE_s*_rho*_Cp*u_phi[i][qp]*grad_T(0)*JxW[qp];
          Fv(i) += -tau_E*RE_s*_rho*_Cp*u_phi[i][qp]*grad_T(1)*JxW[qp];
          if( this->_dim == 3 )
          {
          (*Fw)(i) += -tau_E*RE_s*_rho*_Cp*u_phi[i][qp]*grad_T(2)*JxW[qp];
          }
          }
        */
  
        for (unsigned int i=0; i != n_T_dofs; i++)
          {
            FT(i) += -tau_E*RE_s*( _rho*_Cp*U*T_gradphi[i][qp]
                                  + _k*(T_hessphi[i][qp](0,0) + T_hessphi[i][qp](1,1) + T_hessphi[i][qp](2,2)) 
                                  )*JxW[qp];
            if (compute_jacobian)
              {
                for (unsigned int j=0; j != n_T_dofs; ++j)
                  {
                    KTT(i,j) += -tau_E*
                      (d_RE_s_dT*T_phi[j][qp] +
                       d_RE_s_dgradT*T_gradphi[j][qp] +
                       d_RE_s_dhessT.contract(T_hessphi[j][qp])
                      ) *
                      ( _rho*_Cp*U*T_gradphi[i][qp]
                        + _k*(T_hessphi[i][qp](0,0) +
                              T_hessphi[i][qp](1,1) +
                              T_hessphi[i][qp](2,2)) 
                      )*JxW[qp]
                      * context.get_fixed_solution_derivative();
                  }
                for (unsigned int j=0; j != n_u_dofs; ++j)
                  {
                    KTu(i,j) += -tau_E*RE_s*
                      ( _rho*_Cp*u_phi[j][qp]*T_gradphi[i][qp](0) )*JxW[qp]
                      * context.get_fixed_solution_derivative();
                    KTu(i,j) +=
                      -(tau_E*d_RE_s_dU(0)+d_tau_E_dU(0)*RE_s)*u_phi[j][qp]*
                      ( _rho*_Cp*U*T_gradphi[i][qp]
                        + _k*(T_hessphi[i][qp](0,0) +
                              T_hessphi[i][qp](1,1) +
                              T_hessphi[i][qp](2,2)) 
                      )*JxW[qp]
                      * context.get_fixed_solution_derivative();
                    KTv(i,j) += -tau_E*RE_s*
                      ( _rho*_Cp*u_phi[j][qp]*T_gradphi[i][qp](1) )*JxW[qp]
                      * context.get_fixed_solution_derivative();
                    KTv(i,j) +=
                      -(tau_E*d_RE_s_dU(1)+d_tau_E_dU(1)*RE_s)*u_phi[j][qp]*
                      ( _rho*_Cp*U*T_gradphi[i][qp]
                        + _k*(T_hessphi[i][qp](0,0) +
                              T_hessphi[i][qp](1,1) +
                              T_hessphi[i][qp](2,2)) 
                      )*JxW[qp]
                      * context.get_fixed_solution_derivative();
                  }
                if(this->_dim == 3)
                  {
                    for (unsigned int j=0; j != n_u_dofs; ++j)
                      {
                        (*KTw)(i,j) += -tau_E*RE_s*
                          ( _rho*_Cp*u_phi[j][qp]*T_gradphi[i][qp](2) )*JxW[qp]
                          * context.get_fixed_solution_derivative();
                        (*KTw)(i,j) +=
                          -(tau_E*d_RE_s_dU(2)+d_tau_E_dU(2)*RE_s)*u_phi[j][qp]*
                          ( _rho*_Cp*U*T_gradphi[i][qp]
                            + _k*(T_hessphi[i][qp](0,0) +
                                  T_hessphi[i][qp](1,1) +
                                  T_hessphi[i][qp](2,2)) 
                          )*JxW[qp]
                          * context.get_fixed_solution_derivative();
                      }
                  }
              }
          }
      }

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("HeatTransferAdjointStabilization::element_time_derivative");
#endif
    return;
  }

  void HeatTransferAdjointStabilization::mass_residual( bool compute_jacobian,
                                                        AssemblyContext& context,
                                                        CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("HeatTransferAdjointStabilization::mass_residual");
#endif

    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T_var()).size();
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u_var()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_temp_vars.T_var())->get_JxW();

    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u_var())->get_phi();

    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_temp_vars.T_var())->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(this->_temp_vars.T_var())->get_dphi();

    const std::vector<std::vector<libMesh::RealTensor> >& T_hessphi =
      context.get_element_fe(this->_temp_vars.T_var())->get_d2phi();

    /*
      const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u_var()).size();

      const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u_var())->get_phi();

      libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u_var()); // R_{p}
      libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v_var()); // R_{p}
      libMesh::DenseSubVector<libMesh::Number> *Fw = NULL;
      if(this->_dim == 3)
      {
      Fw = &context.get_elem_residual(this->_flow_vars.w_var()); // R_{w}
      }
    */

    libMesh::DenseSubVector<libMesh::Number> &FT = context.get_elem_residual(this->_temp_vars.T_var()); // R_{T}
    libMesh::DenseSubMatrix<libMesh::Number> &KTT = 
      context.get_elem_jacobian(_temp_vars.T_var(), _temp_vars.T_var()); // J_{TT}
    libMesh::DenseSubMatrix<libMesh::Number> &KTu = 
      context.get_elem_jacobian(_temp_vars.T_var(), _flow_vars.u_var()); // J_{Tu}
    libMesh::DenseSubMatrix<libMesh::Number> &KTv = 
      context.get_elem_jacobian(_temp_vars.T_var(), _flow_vars.v_var()); // J_{Tv}
    libMesh::DenseSubMatrix<libMesh::Number> *KTw = NULL;

    if(this->_dim == 3)
      {
        KTw = &context.get_elem_jacobian
          (_temp_vars.T_var(), _flow_vars.w_var()); // J_{Tw}
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::FEBase* fe = context.get_element_fe(this->_temp_vars.T_var());

        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::RealGradient U( context.fixed_interior_value( this->_flow_vars.u_var(), qp ),
                                 context.fixed_interior_value( this->_flow_vars.v_var(), qp ) );
        if( this->_dim == 3 )
          U(2) = context.fixed_interior_value( this->_flow_vars.w_var(), qp );
      
        //libMesh::RealGradient grad_T = context.fixed_interior_gradient( this->_temp_vars.T_var(), qp );

        libMesh::Real tau_E, RE_t;
        libMesh::Real d_tau_E_d_rho, d_RE_t_dT;
        libMesh::Gradient d_tau_E_dU;
       
        if (compute_jacobian)
          {
            this->_stab_helper.compute_tau_energy_and_derivs
              ( context, G, _rho, _Cp, _k,  U,
                tau_E, d_tau_E_d_rho, d_tau_E_dU, false );
            this->_stab_helper.compute_res_energy_transient_and_derivs
              ( context, qp, _rho, _Cp,
                RE_t, d_RE_t_dT );
          }
        else
          {
            tau_E = this->_stab_helper.compute_tau_energy
              ( context, G, _rho, _Cp, _k,  U, false );

            RE_t = this->_stab_helper.compute_res_energy_transient
              ( context, qp, _rho, _Cp );
          }


        /*
          for (unsigned int i=0; i != n_u_dofs; i++)
          {
          Fu(i) += -tau_E*RE_t*_rho*_Cp*u_phi[i][qp]*grad_T(0)*JxW[qp];
          Fv(i) += -tau_E*RE_t*_rho*_Cp*u_phi[i][qp]*grad_T(1)*JxW[qp];
          if( this->_dim == 3 )
          {
          (*Fw)(i) += -tau_E*RE_t*_rho*_Cp*u_phi[i][qp]*grad_T(2)*JxW[qp];
          }
          }
        */

        for (unsigned int i=0; i != n_T_dofs; i++)
          {
            FT(i) += tau_E*RE_t*( _rho*_Cp*U*T_gradphi[i][qp]
                                  + _k*(T_hessphi[i][qp](0,0) + T_hessphi[i][qp](1,1) + T_hessphi[i][qp](2,2)) 
                                  )*JxW[qp];
            if (compute_jacobian)
              {
                for (unsigned int j=0; j != n_T_dofs; ++j)
                  {
                    KTT(i,j) += 
                      (tau_E*d_RE_t_dT)*T_phi[j][qp]*
                      ( _rho*_Cp*U*T_gradphi[i][qp]
                        + _k*(T_hessphi[i][qp](0,0) +
                              T_hessphi[i][qp](1,1) +
                              T_hessphi[i][qp](2,2)) 
                      )*JxW[qp];
                  }
                for (unsigned int j=0; j != n_u_dofs; ++j)
                  {
                    KTu(i,j) +=
                      d_tau_E_dU(0)*u_phi[j][qp]*RE_t*
                      ( _rho*_Cp*U*T_gradphi[i][qp]
                        + _k*(T_hessphi[i][qp](0,0) +
                              T_hessphi[i][qp](1,1) +
                              T_hessphi[i][qp](2,2)) 
                      )*JxW[qp];
                    KTu(i,j) +=
                      tau_E*RE_t*
                      ( _rho*_Cp*u_phi[j][qp]*T_gradphi[i][qp](0)
                        + _k*(T_hessphi[i][qp](0,0) +
                              T_hessphi[i][qp](1,1) +
                              T_hessphi[i][qp](2,2)) 
                      )*JxW[qp];
                    KTv(i,j) +=
                      d_tau_E_dU(1)*u_phi[j][qp]*RE_t*
                      ( _rho*_Cp*U*T_gradphi[i][qp]
                        + _k*(T_hessphi[i][qp](0,0) +
                              T_hessphi[i][qp](1,1) +
                              T_hessphi[i][qp](2,2)) 
                      )*JxW[qp];
                    KTv(i,j) +=
                      tau_E*RE_t*
                      ( _rho*_Cp*u_phi[j][qp]*T_gradphi[i][qp](1)
                        + _k*(T_hessphi[i][qp](0,0) +
                              T_hessphi[i][qp](1,1) +
                              T_hessphi[i][qp](2,2)) 
                      )*JxW[qp];
                  }
                if(this->_dim == 3)
                  {
                    for (unsigned int j=0; j != n_u_dofs; ++j)
                      {
                        (*KTw)(i,j) +=
                          d_tau_E_dU(2)*u_phi[j][qp]*RE_t*
                          ( _rho*_Cp*U*T_gradphi[i][qp]
                            + _k*(T_hessphi[i][qp](0,0) +
                                  T_hessphi[i][qp](1,1) +
                                  T_hessphi[i][qp](2,2)) 
                          )*JxW[qp];
                        (*KTw)(i,j) +=
                          tau_E*RE_t*
                          ( _rho*_Cp*u_phi[j][qp]*T_gradphi[i][qp](2)
                            + _k*(T_hessphi[i][qp](0,0) +
                                  T_hessphi[i][qp](1,1) +
                                  T_hessphi[i][qp](2,2)) 
                          )*JxW[qp];
                      }
                  }
              }
          }

      }

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("HeatTransferAdjointStabilization::mass_residual");
#endif
    return;
  }

} // namespace GRINS
