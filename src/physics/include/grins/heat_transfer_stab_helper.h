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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef GRINS_HEAT_TRANSFER_STAB_HELPER_H
#define GRINS_HEAT_TRANSFER_STAB_HELPER_H

//GRINS
#include "grins/stab_helper.h"

// libMesh
#include "libmesh/fem_context.h"

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  class HeatTransferStabilizationHelper : public StabilizationHelper
  {
  public:

    HeatTransferStabilizationHelper( const GetPot& input );

    ~HeatTransferStabilizationHelper();

    libMesh::Real compute_tau_energy( libMesh::FEMContext& c,
				      libMesh::RealTensor& G,
				      libMesh::Real rho,
				      libMesh::Real cp,
				      libMesh::Real k,
				      libMesh::Gradient U,
				      bool is_steady ) const;

  protected:

    libMesh::Real _C, _tau_factor;

  }; // class HeatTransferStabilizationHelper

  /* ------------- Inline Functions ---------------*/
  inline
  libMesh::Real HeatTransferStabilizationHelper::compute_tau_energy( libMesh::FEMContext& c,
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
  
}
#endif // GRINS_HEAT_TRANSFER_STAB_HELPER_H
