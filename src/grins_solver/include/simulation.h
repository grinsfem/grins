//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef SIMULATION_H
#define SIMULATION_H

// libMesh
#include "getpot.h"

// GRINS
#include "visualization.h"

namespace GRINS
{
  class Simulation
  {
  public:
    
    Simulation( const GetPot& input,
		GRINS::PhysicsFactory* physics_factory,
		GRINS::SolverFactory* solver_factory );

    ~Simulation();
	
    void run();

    void output_vis();

  private:

    GRINS::MultiphysicsSystem* _system;

    GRINS::Solver* _solver;

    GRINS::Visualization _vis;

  };
}
#endif // SIMULATION_H
