//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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

// libMesh
#include "libmesh/parallel.h"

// SuspendedCable
#include "suspended_cable_solver_factory.h"

int main(int argc, char* argv[])
{
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

// Initialize libMesh library.
      libMesh::LibMeshInit libmesh_init(argc, argv);

      libMesh::out << "Starting GRINS with command:\n";
      for (int i=0; i != argc; ++i)
        libMesh::out << argv[i] << ' ';
      libMesh::out << std::endl;

      GRINS::SimulationBuilder sim_builder;

std::tr1::shared_ptr<GRINS::SuspendedCableSolverFactory> cable_factory( new GRINS::SuspendedCableSolverFactory );

sim_builder.attach_solver_factory( cable_factory );

      GRINS::Simulation grins( libMesh_inputfile,
                               sim_builder,
                               libmesh_init.comm() );

grins.run();

return 0;
}
