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

#ifndef GRINS_HEAT_TRANSFER_STAB_HELPER_H
#define GRINS_HEAT_TRANSFER_STAB_HELPER_H

//GRINS
#include "grins/stab_helper.h"
#include "grins/assembly_context.h"

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  // Forward declarations
  class PrimitiveTempFEVariables;
  class VelocityVariable;
  class PressureFEVariable;

  class HeatTransferStabilizationHelper : public StabilizationHelper
  {
  public:

    HeatTransferStabilizationHelper( const std::string & helper_name,
                                     const GetPot& input );

    ~HeatTransferStabilizationHelper();

    libMesh::Real compute_res_energy_steady( AssemblyContext& context,
                                             unsigned int qp,
                                             const libMesh::Real rho,
                                             const libMesh::Real Cp,
                                             const libMesh::Real k ) const;


    void compute_res_energy_steady_and_derivs( AssemblyContext& context,
                                               unsigned int qp,
                                               const libMesh::Real rho,
                                               const libMesh::Real Cp,
                                               const libMesh::Real k,
                                               libMesh::Real &res,
                                               libMesh::Real &d_res_dT,
                                               libMesh::Gradient &d_res_dgradT,
                                               libMesh::Tensor   &d_res_dhessT,
                                               libMesh::Gradient &d_res_dU
                                               ) const;

    libMesh::Real compute_res_energy_transient( AssemblyContext& context,
                                                unsigned int qp,
                                                const libMesh::Real rho,
                                                const libMesh::Real Cp
                                                ) const;


    void compute_res_energy_transient_and_derivs( AssemblyContext& context,
                                                  unsigned int qp,
                                                  const libMesh::Real rho,
                                                  const libMesh::Real Cp,
                                                  libMesh::Real &res,
                                                  libMesh::Real &d_res_dTdot
                                                  ) const;

    libMesh::Real compute_tau_energy( AssemblyContext& c,
                                      libMesh::RealTensor& G,
                                      libMesh::Real rho,
                                      libMesh::Real cp,
                                      libMesh::Real k,
                                      libMesh::Gradient U,
                                      bool is_steady ) const;

    void compute_tau_energy_and_derivs( AssemblyContext& c,
                                        libMesh::RealTensor& G,
                                        libMesh::Real rho,
                                        libMesh::Real cp,
                                        libMesh::Real k,
                                        libMesh::Gradient U,
                                        libMesh::Real &tau_E,
                                        libMesh::Real &d_tau_E_d_rho,
                                        libMesh::Gradient &d_tau_E_d_U,
                                        bool is_steady ) const;

  protected:

    libMesh::Real _C, _tau_factor;

    const PrimitiveTempFEVariables& _temp_vars;

    const VelocityVariable& _flow_vars;
    const PressureFEVariable& _press_var;

  }; // class HeatTransferStabilizationHelper

  /* ------------- Inline Functions ---------------*/
  inline
  libMesh::Real HeatTransferStabilizationHelper::compute_tau_energy( AssemblyContext& c,
                                                                     libMesh::RealTensor& G,
                                                                     libMesh::Real rho,
                                                                     libMesh::Real cp,
                                                                     libMesh::Real k,
                                                                     libMesh::Gradient U,
                                                                     bool is_steady ) const
  {
    libMesh::Real tau = (rho*cp*U)*(G*(rho*cp*U)) + this->_C*k*k*G.contract(G);

    if(!is_steady)
      tau += (2.0*rho*cp/c.get_deltat_value())*(2.0*rho*cp/c.get_deltat_value());

    return this->_tau_factor/std::sqrt(tau);
  }


  inline
  void HeatTransferStabilizationHelper::compute_tau_energy_and_derivs
  ( AssemblyContext& c,
    libMesh::RealTensor& G,
    libMesh::Real rho,
    libMesh::Real cp,
    libMesh::Real k,
    libMesh::Gradient U,
    libMesh::Real &tau_E,
    libMesh::Real &d_tau_E_d_rho,
    libMesh::Gradient &d_tau_E_d_U,
    bool is_steady ) const
  {
    libMesh::Gradient rhocpU = rho*cp*U;
    libMesh::Gradient GrhocpU = G*rhocpU;
    libMesh::Real rhocpUGrhocpU = rhocpU * GrhocpU;
    libMesh::Real GG = G.contract(G);
    tau_E = (rhocpU)*(GrhocpU) + this->_C*k*k*GG;
    d_tau_E_d_rho = rhocpUGrhocpU*2/rho/cp;
    d_tau_E_d_U = 2*rho*cp*GrhocpU;

    if(!is_steady)
      {
        libMesh::Real two_rhocp_over_dt = 2*rho*cp/c.get_deltat_value();
        tau_E += two_rhocp_over_dt * two_rhocp_over_dt;
        d_tau_E_d_rho += 4*two_rhocp_over_dt/c.get_deltat_value();
      }

    // But what we've computed so far isn't tau; we need
    // tau = _tau_factor/ sqrt(our_tau)

    libMesh::Real root_oldtau = std::sqrt(tau_E);
    libMesh::Real d_tau_d_oldtau = -this->_tau_factor / (tau_E*root_oldtau) / 2;

    d_tau_E_d_rho         = d_tau_d_oldtau * d_tau_E_d_rho;
    d_tau_E_d_U           = d_tau_d_oldtau * d_tau_E_d_U;

    tau_E = this->_tau_factor / root_oldtau;
  }

}
#endif // GRINS_HEAT_TRANSFER_STAB_HELPER_H
