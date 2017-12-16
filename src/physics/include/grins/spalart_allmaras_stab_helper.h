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

#ifndef GRINS_SPALART_ALLMARAS_STAB_HELPER_H
#define GRINS_SPALART_ALLMARAS_STAB_HELPER_H

// GRINS
#include "grins/stab_helper.h"
#include "grins/assembly_context.h"
#include "grins/spalart_allmaras_helper.h"
#include "grins/spalart_allmaras_parameters.h"
#include "grins/parameter_user.h"

//Utils
#include "grins/distance_function.h"

// libMesh foward declarations
class GetPot;

namespace GRINS
{
  // Forward declarations
  class TurbulenceFEVariables;

  class SpalartAllmarasStabilizationHelper : public StabilizationHelper
  {
  public:

    SpalartAllmarasStabilizationHelper( const std::string& helper_name, const GetPot& input );

    ~SpalartAllmarasStabilizationHelper(){};

    libMesh::Real compute_tau_spalart( AssemblyContext& c,
                                       unsigned int qp,
                                       libMesh::RealGradient& g,
                                       libMesh::RealTensor& G,
                                       libMesh::Real rho,
                                       libMesh::Gradient U,
                                       libMesh::Real mu,
                                       bool is_steady ) const;

    void compute_tau_spalart_and_derivs( AssemblyContext& c,
                                         unsigned int qp,
                                         libMesh::RealGradient& g,
                                         libMesh::RealTensor& G,
                                         libMesh::Real rho,
                                         libMesh::Gradient U,
                                         libMesh::Real T,
                                         libMesh::Real &tau_M,
                                         libMesh::Real &d_tau_M_d_rho,
                                         libMesh::Gradient &d_tau_M_d_U,
                                         bool is_steady ) const;

    libMesh::Real compute_tau( AssemblyContext& c,
                               unsigned int qp,
                               libMesh::Real mat_prop_sq,
                               libMesh::RealGradient& g,
                               libMesh::RealTensor& G,
                               libMesh::Real rho,
                               libMesh::Gradient U,
                               bool is_steady ) const;

    void compute_tau_and_derivs( AssemblyContext& c,
                                 unsigned int qp,
                                 libMesh::Real mat_prop_sq,
                                 libMesh::RealGradient& g,
                                 libMesh::RealTensor& G,
                                 libMesh::Real rho,
                                 libMesh::Gradient U,
                                 libMesh::Real& tau,
                                 libMesh::Real& d_tau_d_rho,
                                 libMesh::Gradient& d_tau_d_U,
                                 bool is_steady ) const;


    libMesh::Real compute_res_spalart_steady( AssemblyContext& context,
                                              unsigned int qp,
                                              const libMesh::Real rho,
                                              const libMesh::Real mu,
                                              const libMesh::Real distance_qp,
                                              const bool infinite_distance) const;

    void compute_res_spalart_steady_and_derivs( AssemblyContext& context,
                                                unsigned int qp,
                                                const libMesh::Real rho,
                                                const libMesh::Real mu,
                                                libMesh::Gradient &res_M,
                                                libMesh::Tensor &d_res_M_dgradp,
                                                libMesh::Tensor &d_res_M_dU,
                                                libMesh::Gradient &d_res_Muvw_dgraduvw,
                                                libMesh::Tensor &d_res_Muvw_dhessuvw
                                                ) const;

    libMesh::Real compute_res_spalart_transient( AssemblyContext& context,
                                                 unsigned int qp,
                                                 const libMesh::Real rho ) const;

    void compute_res_spalart_transient_and_derivs( AssemblyContext& context,
                                                   unsigned int qp,
                                                   const libMesh::Real rho,
                                                   libMesh::RealGradient &res_M,
                                                   libMesh::Real &d_res_Muvw_duvw
                                                   ) const;

    // Registers all parameters in this physics and in its property
    // classes
    virtual void register_parameter
    ( const std::string & param_name,
      libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
      const;

  protected:

    libMesh::Real _C, _tau_factor;

    const VelocityVariable& _flow_vars;
    const PressureFEVariable& _press_var;

    const TurbulenceFEVariables& _turbulence_vars;

    SpalartAllmarasHelper _spalart_allmaras_helper;

    SpalartAllmarasParameters _sa_params;

  }; // class SpalartAllmarasStabilizationHelper

  /* ------------- Inline Functions ---------------*/

  inline
  libMesh::Real SpalartAllmarasStabilizationHelper::compute_tau_spalart( AssemblyContext& c,
                                                                         unsigned int qp,
                                                                         libMesh::RealGradient& g,
                                                                         libMesh::RealTensor& G,
                                                                         libMesh::Real rho,
                                                                         libMesh::Gradient U,
                                                                         libMesh::Real mu,
                                                                         bool is_steady ) const
  {
    return this->compute_tau( c, qp, mu*mu, g, G, rho, U, is_steady );
  }

  inline void
  SpalartAllmarasStabilizationHelper::compute_tau_spalart_and_derivs
  ( AssemblyContext& c,
    unsigned int qp,
    libMesh::RealGradient& g,
    libMesh::RealTensor& G,
    libMesh::Real rho,
    libMesh::Gradient U,
    libMesh::Real mu,
    libMesh::Real& tau_M,
    libMesh::Real &d_tau_M_d_rho,
    libMesh::Gradient &d_tau_M_d_U,
    bool is_steady ) const
  {
    this->compute_tau_and_derivs( c, qp, mu*mu, g, G, rho, U, tau_M,
                                  d_tau_M_d_rho, d_tau_M_d_U,
                                  is_steady );
  }


  inline
  libMesh::Real SpalartAllmarasStabilizationHelper::compute_tau( AssemblyContext& c,
                                                                 unsigned int /*qp*/,
                                                                 libMesh::Real mat_prop_sq,
                                                                 libMesh::RealGradient& /*g*/,
                                                                 libMesh::RealTensor& G,
                                                                 libMesh::Real rho,
                                                                 libMesh::Gradient U,
                                                                 bool is_steady ) const
  {
    libMesh::Real tau = (rho*U)*(G*(rho*U)) + this->_C*mat_prop_sq*G.contract(G);

    if(!is_steady)
      tau += (2.0*rho/c.get_deltat_value())*(2.0*rho/c.get_deltat_value());

    return this->_tau_factor/std::sqrt(tau);
  }

  inline
  void
  SpalartAllmarasStabilizationHelper::compute_tau_and_derivs
  ( AssemblyContext& c,
    unsigned int /*qp*/,
    libMesh::Real mat_prop_sq,
    libMesh::RealGradient& /*g*/, // constant
    libMesh::RealTensor& G, // constant
    libMesh::Real rho,
    libMesh::Gradient U,
    libMesh::Real& tau,
    libMesh::Real& d_tau_d_rho,
    libMesh::Gradient& d_tau_d_U,
    bool is_steady ) const
  {
    libMesh::Gradient rhoU = rho*U;
    libMesh::Gradient GrhoU = G*rhoU;
    libMesh::Real rhoUGrhoU = rhoU * GrhoU;
    libMesh::Real GG = G.contract(G);
    tau = rhoUGrhoU + this->_C*mat_prop_sq*GG;
    d_tau_d_rho = rhoUGrhoU*2/rho;
    d_tau_d_U = 2*rho*GrhoU;

    if(!is_steady)
      {
        libMesh::Real two_rho_over_dt = 2*rho/c.get_deltat_value();
        tau += two_rho_over_dt * two_rho_over_dt;
        d_tau_d_rho += 4*two_rho_over_dt/c.get_deltat_value();
      }

    // But what we've computed so far isn't tau; we need
    // tau = _tau_factor/ sqrt(our_tau)

    libMesh::Real root_oldtau = std::sqrt(tau);
    libMesh::Real d_tau_d_oldtau = -this->_tau_factor / (tau*root_oldtau) / 2;

    d_tau_d_rho         = d_tau_d_oldtau * d_tau_d_rho;
    d_tau_d_U           = d_tau_d_oldtau * d_tau_d_U;

    tau = this->_tau_factor / root_oldtau;
  }

}
#endif // GRINS_SPALART_ALLMARAS_STAB_HELPER_H
