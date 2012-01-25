//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "steady_solver.h"

GRINS::SteadySolver::SteadySolver()
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

void GRINS::SteadySolver::solve()
{
  // GRVY timers contained in here (if enabled)
  _system->solve();
  return;
}
