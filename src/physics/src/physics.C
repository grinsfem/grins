//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
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
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "physics.h"

GRINS::Physics::Physics()
  : _local_variable_map_built(false)
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

void GRINS::Physics::set_time_evolving_vars( libMesh::FEMSystem* system )
{
  return;
}

GRINS::VariableMap GRINS::Physics::get_variable_indices_map()
{
  if( !_local_variable_map_built )
    {
      std::cerr << "Error: Must build local variable map before it can be returned."
		<< std::endl;
      libmesh_error(); //TODO: Do we want libmesh_error as our error handler?
    }

  return _var_map;
}

void GRINS::Physics::register_variable_indices( VariableMap& global_map )
{
  return;
}

void GRINS::Physics::attach_bound_func( const unsigned int bc_id, 
					GRINS::BasePointFuncObj* bound_func )
{
  _bound_funcs[bc_id] = bound_func;

  return;
}

#ifdef USE_GRVY_TIMERS
void GRINS::Physics::attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer )
{
  _timer = grvy_timer;
  return;
}
#endif
