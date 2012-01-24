//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "simulation.h"

Simulation::Simulation( const GetPot& input,
			GRINS::PhysicsFactory* physics_factory,
			GRINS::SolverFactory* solver_factory )
  : _solver( solver_factory->build() )
{
  
  return;
}

Simulation::~Simulation()
{
  return;
}

void Simulation::run()
{
  this->print_screen_info();
  
  _solver->solve();

  return;
}
