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
#ifndef INC_NAVIER_STOKES_STAB_BASE_H
#define INC_NAVIER_STOKES_STAB_BASE_H

//GRINS
#include "grins/inc_navier_stokes_base.h"
#include "grins/inc_navier_stokes_stab_helper.h"

//! GRINS namespace
namespace GRINS
{
  class IncompressibleNavierStokesStabilizationBase : public IncompressibleNavierStokesBase
  {

  public:

    IncompressibleNavierStokesStabilizationBase( const GRINS::PhysicsName& physics_name, const GetPot& input );

    virtual ~IncompressibleNavierStokesStabilizationBase();

    //! Initialize context for added physics variables
    virtual void init_context( libMesh::FEMContext& context );

    libMesh::Real compute_res_continuity( libMesh::FEMContext& context,
					  unsigned int qp ) const;
    
    libMesh::RealGradient compute_res_momentum_steady( libMesh::FEMContext& context,
						       unsigned int qp ) const;
    
    libMesh::RealGradient compute_res_momentum_transient( libMesh::FEMContext& context,
							  unsigned int qp ) const;

  protected:

    IncompressibleNavierStokesStabilizationHelper _stab_helper;
    
  private:

    IncompressibleNavierStokesStabilizationBase();

  }; // End IncompressibleNavierStokesStabilizationBase class declarations

} // End namespace GRINS

#endif //INC_NAVIER_STOKES_STAB_BASE_H
