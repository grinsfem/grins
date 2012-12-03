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

#include "grins_steady_solver.h"

namespace GRINS
{

  SteadySolver::SteadySolver( const GetPot& input )
    : Solver( input )
  {
    return;
  }

  SteadySolver::~SteadySolver()
  {
    return;
  }

  void SteadySolver::init_time_solver(MultiphysicsSystem* system)
  {
    libMesh::SteadySolver* time_solver = new libMesh::SteadySolver( *(system) );

    system->time_solver = AutoPtr<TimeSolver>(time_solver);
    return;
  }

  void SteadySolver::solve( MultiphysicsSystem* system,
			    std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
			    std::tr1::shared_ptr<Visualization> vis,
			    bool output_vis, 
			    bool output_residual )
  {
    // GRVY timers contained in here (if enabled)
    system->solve();

    if( output_vis ) vis->output( equation_system );

    if( output_residual ) vis->output_residual( equation_system, system );

    return;
  }

} // namespace GRINS
