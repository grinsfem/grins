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

// GRINS
#include "simulation_builder.h"
#include "simulation.h"
#include "model_adaptive_simulation.h"

// GRVY
#ifdef HAVE_GRVY
#include "grvy.h"
#endif

// libMesh
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
 
  if( libMesh_inputfile( "Adaptivity/model_adaptive_sim", false ) )
  {
    std::string adj_input_filename = libMesh_inputfile( "Adaptivity/adjoint_input_filename",
        libMesh_input_filename ); // adjoint input defaults to forward input
    std::string res_input_filename = libMesh_inputfile( "Adaptivity/residual_input_filename",
        adj_input_filename ); // residual input defaults to adjoint input
    GetPot adjoint_input( adj_input_filename );
    GetPot residual_input( res_input_filename );

    GRINS::SimulationBuilder for_sim_builder;
    GRINS::SimulationBuilder adj_sim_builder;
    GRINS::SimulationBuilder res_sim_builder;

    GRINS::ModelAdaptiveSimulation grins(
        libMesh_inputfile,
        adj_input_filename,
        res_input_filename,
        for_sim_builder,
        adj_sim_builder,
        res_sim_builder );

#ifdef USE_GRVY_TIMERS
    grvy_timer.EndTimer("Initialize Solver");

    // Attach GRVY timer to solver
    grins.attach_grvy_timer( &grvy_timer );
#endif
    grins.run();
  }
  else
  {
    GRINS::SimulationBuilder sim_builder;

    GRINS::Simulation grins( libMesh_inputfile, sim_builder );

#ifdef USE_GRVY_TIMERS
    grvy_timer.EndTimer("Initialize Solver");

    // Attach GRVY timer to solver
    grins.attach_grvy_timer( &grvy_timer );
#endif
    grins.run();
  }

#ifdef USE_GRVY_TIMERS
  grvy_timer.Finalize();
 
  if( Parallel::Communicator_World.rank() == 0 ) grvy_timer.Summarize();
#endif

  return 0;
}
