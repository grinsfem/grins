//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2010 The PECOS Development Team
//
// Please see http://pecos.ices.utexas.edu for more information.
//
// This file is part of GRINS.
//
// GRINS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GRINS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GRINS.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------
//
// Driver code for GRINS (FEM solver).
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "grins_mesh_manager.h"
#include "grins_solver.h"

#include <iostream>
#include "grvy.h"

// System types that we might want to instantiate
#include "low_mach_num_navier_stokes_sys.h"

int main(int argc, char* argv[]) {

  // Check command line count.
  if( argc != 2 )
    {
      std::cout << "Error: Must specify libMesh input file." << std::endl;
      exit(1);
    }

  GRVY::GRVY_Timer_Class grvy_timer;
  grvy_timer.Init("GRVY Timer in GRINS - Generic Driver");

  // Initialize libMesh library.
  LibMeshInit libmesh_init(argc, argv);
  
  // Create mesh manager object.
  GRINS::MeshManager meshmanager;
  
  // Create solver object.
  std::string solver_dummy_options = "TODO: Delete me when agreed on constructor arguments.";
  GRINS::Solver<GRINS::LowMachNumberNavierStokesSystem> solver( solver_dummy_options );

  grvy_timer.BeginTimer("Generic Driver - input reading block timing");
  { // Artificial block to destroy objects associated with reading the input once we've read it in.

    // libMesh input file should be first argument
    std::string libMesh_input = argv[1];
    
    // Create our GetPot object. TODO: Finalize decision of GRVY vs. GetPot input.
    GetPot libMesh_inputfile( libMesh_input );
    
    solver.read_input_options( libMesh_inputfile );

    meshmanager.read_input_options( libMesh_inputfile );
  } //Should be done reading input, so we kill the GetPot object.
  grvy_timer.EndTimer("Generic Driver - input reading block timing");

  // TODO: meshmanager.init, meshmanager.read, etc.

  // pass libMesh::Mesh object from meshmanager to solver
  solver.set_mesh( meshmanager.get_mesh() );

  grvy_timer.Finalize();
  grvy_timer.Summarize();

  return 0;
}
