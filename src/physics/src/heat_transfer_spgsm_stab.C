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
#include "grins/heat_transfer_spgsm_stab.h"

// libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{

  HeatTransferSPGSMStabilization::HeatTransferSPGSMStabilization( const std::string& physics_name, 
                                                                  const GetPot& input )
    : HeatTransferStabilizationBase(physics_name,input)
  {
    return;
  }

  HeatTransferSPGSMStabilization::~HeatTransferSPGSMStabilization()
  {
    return;
  }

  void HeatTransferSPGSMStabilization::element_time_derivative( bool compute_jacobian,
                                                                AssemblyContext& context,
                                                                CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("HeatTransferSPGSMStabilization::element_time_derivative");
#endif

    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T_var()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_temp_vars.T_var())->get_JxW();

    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(this->_temp_vars.T_var())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &FT = context.get_elem_residual(this->_temp_vars.T_var()); // R_{T}

    libMesh::FEBase* fe = context.get_element_fe(this->_temp_vars.T_var());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::RealGradient U( context.interior_value( this->_flow_vars.u_var(), qp ),
                                 context.interior_value( this->_flow_vars.v_var(), qp ) );
        if( this->_dim == 3 )
          {
            U(2) = context.interior_value( this->_flow_vars.w_var(), qp );
          }
      
        libMesh::Real tau_E = this->_stab_helper.compute_tau_energy( context, G, _rho, _Cp, _k,  U, this->_is_steady );

        libMesh::Real RE_s = this->_stab_helper.compute_res_energy_steady( context, qp, _rho, _Cp, _k );

        for (unsigned int i=0; i != n_T_dofs; i++)
          {
            FT(i) += -tau_E*RE_s*_rho*_Cp*U*T_gradphi[i][qp]*JxW[qp];
          }

        if( compute_jacobian )
          {
            libmesh_not_implemented();
          }

      }

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("HeatTransferSPGSMStabilization::element_time_derivative");
#endif
    return;
  }

  void HeatTransferSPGSMStabilization::mass_residual( bool /*compute_jacobian*/,
                                                      AssemblyContext& context,
                                                      CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("HeatTransferSPGSMStabilization::mass_residual");
#endif

    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T_var()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_temp_vars.T_var())->get_JxW();

    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(this->_temp_vars.T_var())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &FT = context.get_elem_residual(this->_temp_vars.T_var()); // R_{T}

    libMesh::FEBase* fe = context.get_element_fe(this->_temp_vars.T_var());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::RealGradient U( context.fixed_interior_value( this->_flow_vars.u_var(), qp ),
                                 context.fixed_interior_value( this->_flow_vars.v_var(), qp ) );
        if( this->_dim == 3 )
          {
            U(2) = context.fixed_interior_value( this->_flow_vars.w_var(), qp );
          }
      
        libMesh::Real tau_E = this->_stab_helper.compute_tau_energy( context, G, _rho, _Cp, _k,  U, false );

        libMesh::Real RE_t = this->_stab_helper.compute_res_energy_transient( context, qp, _rho, _Cp );

        for (unsigned int i=0; i != n_T_dofs; i++)
          {
            FT(i) += tau_E*RE_t*_rho*_Cp*U*T_gradphi[i][qp]*JxW[qp];
          }

      }

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("HeatTransferSPGSMStabilization::mass_residual");
#endif
    return;
  }

} // namespace GRINS
