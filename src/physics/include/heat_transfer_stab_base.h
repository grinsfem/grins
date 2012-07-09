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
#ifndef HEAT_TRANSFER_STAB_BASE_H
#define HEAT_TRANSFER_STAB_BASE_H

//GRINS
#include "heat_transfer_base.h"
#include "heat_transfer_stab_helper.h"

//! GRINS namespace
namespace GRINS
{
  class HeatTransferStabilizationBase : public HeatTransferBase
  {

  public:

    HeatTransferStabilizationBase( const GRINS::PhysicsName& physics_name, const GetPot& input );

    virtual ~HeatTransferStabilizationBase();

    //! Read options from GetPot input file. By default, nothing is read.
    virtual void read_input_options( const GetPot& input );

    //! Initialize context for added physics variables
    virtual void init_context( libMesh::DiffContext &context );
    
    libMesh::Real compute_res_steady( libMesh::FEMContext& context,
				      unsigned int qp ) const;
    
    libMesh::Real compute_res_transient( libMesh::FEMContext& context,
					 unsigned int qp ) const;

  protected:

    HeatTransferStabilizationHelper _stab_helper;
    
  private:

    HeatTransferStabilizationBase();

  }; // End HeatTransferStabilizationBase class declarations

} // End namespace GRINS

#endif //HEAT_TRANSFER_STAB_BASE_H
