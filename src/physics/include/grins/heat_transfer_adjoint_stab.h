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
#ifndef GRINS_HEAT_TRANSFER_ADJOINT_STAB_H
#define GRINS_HEAT_TRANSFER_ADJOINT_STAB_H

//GRINS
#include "grins/heat_transfer_stab_base.h"

//! GRINS namespace
namespace GRINS
{
  //! Adds VMS-based stabilization to LowMachNavierStokes physics class
  class HeatTransferAdjointStabilization : public HeatTransferStabilizationBase
  {

  public:

    HeatTransferAdjointStabilization( const GRINS::PhysicsName& physics_name, const GetPot& input );
    virtual ~HeatTransferAdjointStabilization();

    virtual void element_time_derivative( bool compute_jacobian,
					  libMesh::FEMContext& context,
					  CachedValues& cache );

    virtual void mass_residual( bool compute_jacobian,
				libMesh::FEMContext& context,
				CachedValues& cache );
    
  private:
    HeatTransferAdjointStabilization();

  }; // End HeatTransferAdjointStabilization class declarations

} // End namespace GRINS

#endif // GRINS_HEAT_TRANSFER_ADJOINT_STAB_H
