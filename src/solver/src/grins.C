//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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


#include "grins_config.h"

#include <iostream>

// GRINS
#include "grins/simulation_builder.h"
#include "grins/simulation.h"

// GRVY
#ifdef GRINS_HAVE_GRVY
#include "grvy.h"
#endif

// libMesh
#include "libmesh/parallel.h"

int main(int argc, char* argv[])
{
  /* Echo GRINS version, libMesh version, and command */
  libMesh::out << "=========================================================="
               << std::endl;
  libMesh::out << "GRINS Version: " << GRINS_BUILD_VERSION << std::endl
               << "libMesh Version: " << LIBMESH_BUILD_VERSION << std::endl
               << "Running with command:\n";

  for (int i=0; i != argc; ++i)
    std::cout << argv[i] << ' ';

  std::cout << std::endl
            << "=========================================================="
            << std::endl;

#ifdef GRINS_USE_GRVY_TIMERS
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
  
  // Initialize libMesh library.
  libMesh::LibMeshInit libmesh_init(argc, argv);

  // Create our GetPot object.
  GetPot libMesh_inputfile( libMesh_input_filename );

  GetPot command_line(argc,argv);

  // GetPot doesn't throw an error for a nonexistent file?
  {
    std::ifstream i(libMesh_input_filename.c_str());
    if (!i)
      {
        std::cerr << "Error: Could not read from libMesh input file "
                << libMesh_input_filename << std::endl;
        exit(1);
      }
  }

#ifdef GRINS_USE_GRVY_TIMERS
  grvy_timer.BeginTimer("Initialize Solver");
#endif

  GRINS::SimulationBuilder sim_builder;

  GRINS::Simulation grins( libMesh_inputfile,
                           command_line,
			   sim_builder,
                           libmesh_init.comm() );

#ifdef GRINS_USE_GRVY_TIMERS
  grvy_timer.EndTimer("Initialize Solver");

  // Attach GRVY timer to solver
  grins.attach_grvy_timer( &grvy_timer );
#endif

grins.run();

#ifdef GRINS_USE_GRVY_TIMERS
  grvy_timer.Finalize();
 
  if( Parallel::Communicator_World.rank() == 0 ) grvy_timer.Summarize();
#endif

  return 0;
}
