//-----------------------------------------------------------------------bl-
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

void GRINS::SteadySolver::solve( GRINS::Visualization* vis )
{
  // GRVY timers contained in here (if enabled)
  _system->solve();

  // Output solution and residual, if requested at runtime.
  vis->output();
  vis->output_residual();

  return;
}
