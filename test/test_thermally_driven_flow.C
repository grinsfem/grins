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
#include "grins_mesh_manager.h"
#include "grins_solver.h"

// System types that we might want to instantiate
#include "multiphysics_sys.h"

#include "exact_solution.h"

int main(int argc, char* argv[]) 
{

  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify libMesh input file." << std::endl;
      exit(1); // TODO: something more sophisticated for parallel runs?
    }

  // Initialize libMesh library.
  LibMeshInit libmesh_init(argc, argv);
  
  // Create mesh manager object.
  GRINS::MeshManager meshmanager;
  
  // Create solver object.
  GRINS::Solver<GRINS::MultiphysicsSystem> solver;

  // Filename of file where comparison solution is stashed
  std::string solution_file;

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

    solution_file = libMesh_inputfile( "ExactSolution/solution_file", "DIE!" );

  } //Should be done reading input, so we kill the GetPot object.

  // Do solve here
  solver.solve();

  // Get equation systems to create ExactSolution object
  GRINS::MultiphysicsSystem* system = solver.get_system();

  EquationSystems & es = system->get_equation_systems ();

  // Create Exact solution object and attach exact solution quantities
  ExactSolution exact_sol(es);
  
  EquationSystems es_ref( *meshmanager.get_mesh() );

  es_ref.read( solution_file );

  exact_sol.attach_reference_solution( &es_ref );
  
  // Compute error and get it in various norms
  exact_sol.compute_error("GRINS", "u");
  exact_sol.compute_error("GRINS", "v");

  if( (meshmanager.get_mesh())->mesh_dimension() == 3 )
    exact_sol.compute_error("GRINS", "w");

  exact_sol.compute_error("GRINS", "p");
  exact_sol.compute_error("GRINS", "T");

  double u_l2error = exact_sol.l2_error("GRINS", "u");
  double u_h1error = exact_sol.h1_error("GRINS", "u");

  double v_l2error = exact_sol.l2_error("GRINS", "v");
  double v_h1error = exact_sol.h1_error("GRINS", "v");

  double p_l2error = exact_sol.l2_error("GRINS", "p");
  double p_h1error = exact_sol.h1_error("GRINS", "p");

  double T_l2error = exact_sol.l2_error("GRINS", "T");
  double T_h1error = exact_sol.h1_error("GRINS", "T");
  
  double w_l2error = 0.0, 
         w_h1error = 0.0;

  if( (meshmanager.get_mesh())->mesh_dimension() == 3 )
    {
      w_l2error = exact_sol.l2_error("GRINS", "w");
      w_h1error = exact_sol.h1_error("GRINS", "w");
    }

  int return_flag = 0;

  // This is the tolerance of the iterative linear solver so
  // it's unreasonable to expect anything better than this.
  double tol = 6.0e-12;
  
  if( u_l2error > tol || u_h1error > tol ||
      v_l2error > tol || v_h1error > tol ||
      w_l2error > tol || w_h1error > tol ||
      p_l2error > tol || p_h1error > tol ||
      T_l2error > tol || T_h1error > tol   )
    {
      return_flag = 1;

      std::cout << "Tolerance exceeded for thermally driven flow test." << std::endl
		<< "tolerance = " << tol << std::endl
		<< "u l2 error = " << u_l2error << std::endl
		<< "u h1 error = " << u_h1error << std::endl
		<< "v l2 error = " << v_l2error << std::endl
		<< "v h1 error = " << v_h1error << std::endl
		<< "w l2 error = " << w_l2error << std::endl
		<< "w h1 error = " << w_h1error << std::endl
		<< "p l2 error = " << p_l2error << std::endl
		<< "p h1 error = " << p_h1error << std::endl
		<< "T l2 error = " << T_l2error << std::endl
		<< "T h1 error = " << T_h1error << std::endl;
    }

 return return_flag;
}
