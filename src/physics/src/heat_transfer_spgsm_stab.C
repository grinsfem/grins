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
#include "grins/assembly_context.h"
#include "grins/heat_transfer_macros.h"
#include "grins/heat_transfer_spgsm_stab.h"

// libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{

  template<class K>
  HeatTransferSPGSMStabilization<K>::HeatTransferSPGSMStabilization( const std::string& physics_name,
                                                                     const GetPot& input )
    : HeatTransferStabilizationBase<K>(physics_name,input)
  {
    return;
  }

  template<class K>
  HeatTransferSPGSMStabilization<K>::~HeatTransferSPGSMStabilization()
  {
    return;
  }

  template<class K>
  void HeatTransferSPGSMStabilization<K>::element_time_derivative
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_temp_vars.T())->get_JxW();

    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(this->_temp_vars.T())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &FT = context.get_elem_residual(this->_temp_vars.T()); // R_{T}

    libMesh::FEBase* fe = context.get_element_fe(this->_temp_vars.T());

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

        // Compute Conductivity at this qp
        libMesh::Real _k_qp = this->_k(context, qp);

        libMesh::Real tau_E = this->_stab_helper.compute_tau_energy( context, G, this->_rho, this->_Cp, _k_qp,  U, this->_is_steady );

        libMesh::Real RE_s = this->_stab_helper.compute_res_energy_steady( context, qp, this->_rho, this->_Cp, _k_qp );

        for (unsigned int i=0; i != n_T_dofs; i++)
          {
            FT(i) += -tau_E*RE_s*this->_rho*this->_Cp*U*T_gradphi[i][qp]*JxW[qp];
          }

        if( compute_jacobian )
          {
            libmesh_not_implemented();
          }

      }
  }

  template<class K>
  void HeatTransferSPGSMStabilization<K>::mass_residual( bool compute_jacobian,
                                                         AssemblyContext & context )
  {
    if( compute_jacobian )
      libmesh_not_implemented();

    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_temp_vars.T())->get_JxW();

    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(this->_temp_vars.T())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &FT = context.get_elem_residual(this->_temp_vars.T()); // R_{T}

    libMesh::FEBase* fe = context.get_element_fe(this->_temp_vars.T());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::RealGradient U( context.fixed_interior_value( this->_flow_vars.u(), qp ),
                                 context.fixed_interior_value( this->_flow_vars.v(), qp ) );
        if( this->_flow_vars.dim() == 3 )
          {
            U(2) = context.fixed_interior_value( this->_flow_vars.w(), qp );
          }

        // Compute Conductivity at this qp
        libMesh::Real _k_qp = this->_k(context, qp);

        libMesh::Real tau_E = this->_stab_helper.compute_tau_energy( context, G, this->_rho, this->_Cp, _k_qp,  U, false );

        libMesh::Real RE_t = this->_stab_helper.compute_res_energy_transient( context, qp, this->_rho, this->_Cp );

        for (unsigned int i=0; i != n_T_dofs; i++)
          {
            FT(i) -= tau_E*RE_t*this->_rho*this->_Cp*U*T_gradphi[i][qp]*JxW[qp];
          }

      }
  }

} // namespace GRINS

// Instantiate
INSTANTIATE_HEAT_TRANSFER_SUBCLASS(HeatTransferSPGSMStabilization);
