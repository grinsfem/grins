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

#ifndef GRINS_INC_NAVIER_STOKES_STAB_HELPER_H
#define GRINS_INC_NAVIER_STOKES_STAB_HELPER_H

// GRINS
#include "grins/stab_helper.h"
#include "grins/assembly_context.h"
#include "grins/primitive_flow_variables.h"

// libMesh foward declarations
class GetPot;

namespace GRINS
{
  class IncompressibleNavierStokesStabilizationHelper : public StabilizationHelper
  {
  public:

    IncompressibleNavierStokesStabilizationHelper( const GetPot& input );

    ~IncompressibleNavierStokesStabilizationHelper();

    void init( libMesh::FEMSystem& system );

    libMesh::Real compute_tau_continuity( libMesh::Real tau_M,
                                          libMesh::RealGradient& g  ) const;

    libMesh::Real compute_tau_momentum( AssemblyContext& c,
                                        unsigned int qp,
                                        libMesh::RealGradient& g,
                                        libMesh::RealTensor& G,
                                        libMesh::Real rho,
                                        libMesh::Gradient U,
                                        libMesh::Real T,
                                        bool is_steady ) const;

    libMesh::Real compute_tau( AssemblyContext& c,
                               unsigned int qp,
                               libMesh::Real mat_prop_sq,
                               libMesh::RealGradient& g,
                               libMesh::RealTensor& G,
                               libMesh::Real rho,
                               libMesh::Gradient U,
                               bool is_steady ) const;

    libMesh::Real compute_res_continuity( AssemblyContext& context,
                                          unsigned int qp ) const;

    libMesh::RealGradient compute_res_momentum_steady( AssemblyContext& context,
                                                       unsigned int qp,
                                                       const libMesh::Real rho,
                                                       const libMesh::Real mu ) const;

    libMesh::RealGradient compute_res_momentum_transient( AssemblyContext& context,
                                                          unsigned int qp,
                                                          const libMesh::Real rho ) const;

    /*! \todo Should we inline this? */
    libMesh::RealGradient UdotGradU( libMesh::Gradient& U, libMesh::Gradient& grad_u, 
                                     libMesh::Gradient& grad_v ) const;
    
    /*! \todo Should we inline this? */
    libMesh::RealGradient UdotGradU( libMesh::Gradient& U, libMesh::Gradient& grad_u, 
                                     libMesh::Gradient& grad_v, libMesh::Gradient& grad_w ) const;

    /*! \todo Should we inline this? */
    libMesh::RealGradient div_GradU( libMesh::RealTensor& hess_u, libMesh::RealTensor& hess_v ) const;

    /*! \todo Should we inline this? */
    libMesh::RealGradient div_GradU( libMesh::RealTensor& hess_u, libMesh::RealTensor& hess_v,
                                     libMesh::RealTensor& hess_w ) const;

    /*! \todo Should we inline this? */
    libMesh::RealGradient div_GradU_T( libMesh::RealTensor& hess_u, libMesh::RealTensor& hess_v ) const;

    /*! \todo Should we inline this? */
    libMesh::RealGradient div_GradU_T( libMesh::RealTensor& hess_u, libMesh::RealTensor& hess_v,
                                       libMesh::RealTensor& hess_w ) const;
    
    /*! \todo Should we inline this? */
    libMesh::RealGradient div_divU_I( libMesh::RealTensor& hess_u, libMesh::RealTensor& hess_v ) const;

    /*! \todo Should we inline this? */
    libMesh::RealGradient div_divU_I( libMesh::RealTensor& hess_u, libMesh::RealTensor& hess_v,
                                      libMesh::RealTensor& hess_w ) const;

  protected:

    libMesh::Real _C, _tau_factor;

    PrimitiveFlowVariables _flow_vars;

  }; // class IncompressibleNavierStokesStabilizationHelper

  /* ------------- Inline Functions ---------------*/

  inline
  libMesh::Real IncompressibleNavierStokesStabilizationHelper::compute_tau_continuity( libMesh::Real tau_M,
                                                                                       libMesh::RealGradient& g  ) const
  {
    return this->_tau_factor/(tau_M*(g*g));
  }

  inline
  libMesh::Real IncompressibleNavierStokesStabilizationHelper::compute_tau_momentum( AssemblyContext& c,
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
  
  inline
  libMesh::Real IncompressibleNavierStokesStabilizationHelper::compute_tau( AssemblyContext& c,
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
  IncompressibleNavierStokesStabilizationHelper::compute_tau_and_derivs
    ( AssemblyContext& c,
      unsigned int /*qp*/,
      libMesh::Real mat_prop_sq,
      libMesh::RealGradient& /*g*/, // constant
      libMesh::RealTensor& G, // constant
      libMesh::Real rho,
      libMesh::Gradient U,
      libMesh::Real& tau,
      libMesh::Real& d_tau_d_mat_prop_sq,
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
    d_tau_d_mat_prop_sq = this->_C*GG;
    
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

    d_tau_d_mat_prop_sq = d_tau_d_oldtau * d_tau_d_mat_prop_sq;
    d_tau_d_rho         = d_tau_d_oldtau * d_tau_d_rho;
    d_tau_d_U           = d_tau_d_oldtau * d_tau_d_U;

    tau = this->_tau_factor / roottau;
  }
  
}
#endif // GRINS_INC_NAVIER_STOKES_STAB_HELPER_H
