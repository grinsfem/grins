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

#include "steady_solver.h"

GRINS::SteadySolver::SteadySolver( const GetPot& input )
  : Solver( input )
{
  return;
}

GRINS::SteadySolver::~SteadySolver()
{
  return;
}

void GRINS::SteadySolver::init_time_solver()
{
  libMesh::SteadySolver* time_solver = new libMesh::SteadySolver( *(_system) );

  _system->time_solver = AutoPtr<TimeSolver>(time_solver);
  return;
}

void GRINS::SteadySolver::solve( GRINS::Visualization* vis,
				 bool output_vis, 
				 bool output_residual )
{
  // GRVY timers contained in here (if enabled)
  _system->solve();

  if( output_vis ) vis->output();

  if( output_residual ) vis->output_residual();

  return;
}
