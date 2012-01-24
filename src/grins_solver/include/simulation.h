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
#include "auto_ptr.h"

// GRINS
#include "physics_factory.h"
#include "solver_factory.h"
#include "grins_solver.h"
#include "visualization_factory.h"
#include "visualization.h"

namespace GRINS
{
  class Simulation
  {
  public:
    
    Simulation( const GetPot& input,
		GRINS::PhysicsFactory* physics_factory,
		GRINS::MeshBuilder* mesh_builder,
		GRINS::SolverFactory* solver_factory,
		GRINS::VisualizationFactory* vis_factory );

    ~Simulation();
	
    void run();

    void output_vis();

    void read_input_options( const GetPot& input );

    void print_sim_info();

  private:

    libMesh::AutoPtr<libMesh::Mesh> _mesh;

    libMesh::AutoPtr<libMesh::EquationSystems> _equation_system;

    libMesh::AutoPtr<GRINS::Solver> _solver;

    libMesh::AutoPtr<GRINS::Visualization> _vis;

    // Screen display options
    bool _print_mesh_info;
    bool _print_log_info;
    bool _print_equation_system_info;

  };
}
#endif // SIMULATION_H
