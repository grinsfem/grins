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

#include "physics.h"

GRINS::LowMachNavierStokesVMSStabilization::LowMachNavierStokesVMSStabilization( const std::string& physics_name )
  : Physics(physics_name)
{
  return;
}

GRINS::LowMachNavierStokesVMSStabilization::~LowMachNavierStokesVMSStabilization()
{
  return;
}

void GRINS::LowMachNavierStokesVMSStabilization::read_input_options( const GetPot& input )
{
  return;
}

bool GRINS::LowMachNavierStokesVMSStabilization::element_time_derivative( bool request_jacobian,
									  libMesh::DiffContext& context,
									  libMesh::FEMSystem* system );

