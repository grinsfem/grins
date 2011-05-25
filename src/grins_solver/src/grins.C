//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - a low Mach number Navier-Stokes Finite-Element Solver
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
// Driver code for GRINS (FEM solver).
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "config.h"

#include <iostream>

// GRINS stuff
#include "grins_mesh_manager.h"
#include "grins_solver.h"

#ifdef HAVE_GRVY
// GRVY includes
#include "grvy.h"
#endif

// System types that we might want to instantiate
#include "multiphysics_sys.h"

#include "parallel.h"

int main(int argc, char* argv[]) {

  // Check command line count.
  if( argc != 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify libMesh input file." << std::endl;
      exit(1); // TODO: something more sophisticated for parallel runs?
    }

#ifdef USE_GRVY_TIMERS
  GRVY::GRVY_Timer_Class grvy_timer;
  grvy_timer.Init("GRINS Timer");
#endif

  // Initialize libMesh library.
  LibMeshInit libmesh_init(argc, argv);
  
  // Create mesh manager object.
  GRINS::MeshManager meshmanager;
  
  // Create solver object.
  GRINS::Solver<GRINS::MultiphysicsSystem> solver;

  // Variables we'll want to read in.
  bool output_vis_flag;

#ifdef USE_GRVY_TIMERS
  grvy_timer.BeginTimer("Initialize Solver");
#endif

  { // Artificial block to destroy objects associated with reading the input once we've read it in.

    // libMesh input file should be first argument
    std::string libMesh_input_filename = argv[1];
    
    // Create our GetPot object. TODO: Finalize decision of GRVY vs. GetPot input.
    GetPot libMesh_inputfile( libMesh_input_filename );
    
    // Read solver options
    solver.read_input_options( libMesh_inputfile );

    // Read mesh options
    meshmanager.read_input_options( libMesh_inputfile );

    // Setup and initialize system so system can read it's relavent options
    meshmanager.build_mesh();
    solver.set_mesh( meshmanager.get_mesh() );
    solver.initialize_system( "GRINS", libMesh_inputfile );

    // Read local options
    output_vis_flag = libMesh_inputfile( "vis-options/output_vis_flag", false );

  } //Should be done reading input, so we kill the GetPot object.

#ifdef USE_GRVY_TIMERS
  grvy_timer.EndTimer("Initialize Solver");

  // Attach GRVY timer to solver
  solver.attach_grvy_timer( &grvy_timer );
#endif

  // Do solve here
  solver.solve();

  // Do visualization if we want it.
  if(output_vis_flag) solver.output_visualization(); //TODO: move this in GRINS::Solver

#ifdef USE_GRVY_TIMERS
  grvy_timer.Finalize();
 
  if( Parallel::Communicator_World.rank() == 0 ) grvy_timer.Summarize();
#endif

  return 0;
}
