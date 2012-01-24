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

#include "config.h"

#include <iostream>

// GRINS stuff
#include "mesh_builder.h"
#include "simulation.h"

#ifdef HAVE_GRVY
// GRVY includes
#include "grvy.h"
#endif

#include "parallel.h"

int main(int argc, char* argv[])
{
#ifdef USE_GRVY_TIMERS
  GRVY::GRVY_Timer_Class grvy_timer;
  grvy_timer.Init("GRINS Timer");
#endif

  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify libMesh input file." << std::endl;
      exit(1); // TODO: something more sophisticated for parallel runs?
    }

  // libMesh input file should be first argument
  std::string libMesh_input_filename = argv[1];
  
  // Create our GetPot object.
  GetPot libMesh_inputfile( libMesh_input_filename );

#ifdef USE_GRVY_TIMERS
  grvy_timer.BeginTimer("Initialize Solver");
#endif

  // Initialize libMesh library.
  LibMeshInit libmesh_init(argc, argv);
 
  // MeshBuilder for handling mesh construction
  GRINS::MeshBuilder mesh_builder( libMesh_inputfile );

  // PhysicsFactory handles which GRINS::Physics objects to create
  GRINS::PhysicsFactory physics_factory( libMesh_inputfile );

  // PhysicsFactory handles which GRINS::Solver to use to solve the problem
  GRINS::SolverFactory solver_factory( libMesh_inputfile );

  // VisualizationFactory handles the type of visualization for the simulation
  GRINS::VisualizationFactory vis_factory( );

  GRINS::Simulation grins( libMesh_inputfile,
			   &physics_factory,
			   &mesh_builder,
			   &solver_factory,
			   &vis_factory );

  grins.run();

  grins.output_vis();

#ifdef USE_GRVY_TIMERS
  grvy_timer.EndTimer("Initialize Solver");

  // Attach GRVY timer to solver
  solver.attach_grvy_timer( &grvy_timer );
#endif

#ifdef USE_GRVY_TIMERS
  grvy_timer.Finalize();
 
  if( Parallel::Communicator_World.rank() == 0 ) grvy_timer.Summarize();
#endif

  return 0;
}
