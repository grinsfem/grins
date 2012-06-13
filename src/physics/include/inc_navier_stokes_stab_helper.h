//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef INC_NAVIER_STOKES_STAB_HELPER_H
#define INC_NAVIER_STOKES_STAB_HELPER_H

//libMesh
#include "getpot.h"

//GRINS
#include "stab_helper.h"

namespace GRINS
{
  class IncompressibleNavierStokesStabilizationHelper : public StabilizationHelper
  {
  public:

    IncompressibleNavierStokesStabilizationHelper( const GetPot& input );

    ~IncompressibleNavierStokesStabilizationHelper();

    libMesh::Real compute_tau_continuity( libMesh::Real tau_M,
					  libMesh::RealGradient& g  ) const;

    libMesh::Real compute_tau_momentum( libMesh::FEMContext& c,
					unsigned int qp,
					libMesh::RealGradient& g,
					libMesh::RealTensor& G,
					libMesh::Real rho,
					libMesh::Gradient U,
					libMesh::Real T,
					bool is_steady ) const;

    libMesh::Real compute_tau( libMesh::FEMContext& c,
			       unsigned int qp,
			       libMesh::Real mat_prop_sq,
			       libMesh::RealGradient& g,
			       libMesh::RealTensor& G,
			       libMesh::Real rho,
			       libMesh::Gradient U,
			       bool is_steady ) const;

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

  }; // class IncompressibleNavierStokesStabilizationHelper

  /* ------------- Inline Functions ---------------*/

  inline
  libMesh::Real IncompressibleNavierStokesStabilizationHelper::compute_tau_continuity( libMesh::Real tau_M,
										       libMesh::RealGradient& g  ) const
  {
    return this->_tau_factor/(tau_M*(g*g));
  }

  inline
  libMesh::Real IncompressibleNavierStokesStabilizationHelper::compute_tau_momentum( libMesh::FEMContext& c,
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
  libMesh::Real IncompressibleNavierStokesStabilizationHelper::compute_tau( libMesh::FEMContext& c,
									    unsigned int qp,
									    libMesh::Real mat_prop_sq,
									    libMesh::RealGradient& g,
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
  
}
#endif // INC_NAVIER_STOKES_STAB_HELPER_H
