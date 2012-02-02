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
{
  return;
}

GRINS::Physics::~Physics()
{
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

void GRINS::Physics::attach_dirichlet_bound_func( GRINS::DBCContainer& dirichlet_bcs )
{
  _dirichlet_bound_funcs = dirichlet_bcs;

  return;
}

void GRINS::Physics::attach_neumann_bound_func( GRINS::NBCContainer& neumann_bcs )
{
  _neumann_bound_funcs = neumann_bcs;

  return;
}

#ifdef USE_GRVY_TIMERS
void GRINS::Physics::attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer )
{
  _timer = grvy_timer;
  return;
}
#endif
