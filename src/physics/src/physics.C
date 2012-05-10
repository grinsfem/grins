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

GRINS::Physics::Physics( const std::string& physics_name )
  : _physics_name( physics_name ),
    _bc_handler(NULL)
{
  return;
}

GRINS::Physics::~Physics()
{
  // If a derived class created a bc_handler object, we kill it here.
  if( _bc_handler ) delete _bc_handler;
  return;
}

void GRINS::Physics::read_input_options( const GetPot& input )
{
  return;
}

void GRINS::Physics::set_time_evolving_vars( libMesh::FEMSystem* system )
{
  return;
}

void GRINS::Physics::init_dirichlet_bcs( libMesh::FEMSystem* system )
{
  // Only need to init BC's if the physics actually created a handler
  if( _bc_handler )
    {
      _bc_handler->init_dirichlet_bcs( system );
      _bc_handler->init_dirichlet_bc_func_objs( system );
    }

  return;
}

void GRINS::Physics::attach_neumann_bound_func( GRINS::NBCContainer& neumann_bcs )
{
  _bc_handler->attach_neumann_bound_func( neumann_bcs );
  return;
}

void GRINS::Physics::attach_dirichlet_bound_func( const GRINS::DBCContainer& dirichlet_bc )
{
  _bc_handler->attach_dirichlet_bound_func( dirichlet_bc );
  return;
}

#ifdef USE_GRVY_TIMERS
void GRINS::Physics::attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer )
{
  _timer = grvy_timer;
  return;
}
#endif
