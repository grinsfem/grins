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

// This class
#include "grins/grins_steady_solver.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/solver_context.h"

// libMesh
#include "libmesh/auto_ptr.h"
#include "libmesh/getpot.h"
#include "libmesh/steady_solver.h"


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

  void SteadySolver::solve( SolverContext& context )
  {
    libmesh_assert( context.system );

    if( context.output_vis ) 
      {
	context.postprocessing->update_quantities( *(context.equation_system) );
	context.vis->output( context.equation_system );
      }

    // GRVY timers contained in here (if enabled)
    context.system->solve();

    if( context.output_vis ) 
      {
	context.postprocessing->update_quantities( *(context.equation_system) );
	context.vis->output( context.equation_system );
      }

    if( context.output_residual ) context.vis->output_residual( context.equation_system, context.system );

    return;
  }

} // namespace GRINS
