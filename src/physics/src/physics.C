//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - a low Mach number Navier-Stokes Finite-Element Solver
//
// Copyright (C) 2010,2011 The PECOS Development Team
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
// Definitions for the Physics abstract base class.
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "physics.h"

GRINS::Physics::Physics()
{
  return;
}

GRINS::Physics::~Physics()
{
  return;
}

void GRINS::Physics::read_input_options( GetPot& input )
{
  return;
}

void GRINS::Physics::set_time_evolving_vars( FEMSystem* system )
{
  return;
}

std::map<std::string,GRINS::VariableIndex> GRINS::Physics::get_variable_indices()
{
  return _var_map;
}
