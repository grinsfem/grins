//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
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
#include "mesh_builder.h"
#include "simulation.h"
#include "multiphysics_sys.h"

//libMesh
#include "exact_solution.h"

// GRVY
#ifdef HAVE_GRVY
#include "grvy.h"
#endif

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
  GRINS::PhysicsFactory physics_factory;

  // PhysicsFactory handles which GRINS::Solver to use to solve the problem
  GRINS::SolverFactory solver_factory( libMesh_inputfile );

  // VisualizationFactory handles the type of visualization for the simulation
  GRINS::VisualizationFactory vis_factory( libMesh_inputfile );

  GRINS::Simulation grins( libMesh_inputfile,
			   &physics_factory,
			   &mesh_builder,
			   &solver_factory,
			   &vis_factory );

#ifdef USE_GRVY_TIMERS
  grvy_timer.EndTimer("Initialize Solver");

  // Attach GRVY timer to solver
  grins.attach_grvy_timer( &grvy_timer );
#endif

  // Do solve here
  grins.run();

  // Get equation systems to create ExactSolution object
  std::tr1::shared_ptr<EquationSystems> es = grins.get_equation_system();

  //es->write("foobar.xdr");

  // Create Exact solution object and attach exact solution quantities
  ExactSolution exact_sol(*es);
  
  EquationSystems es_ref( es->get_mesh() );

  // Filename of file where comparison solution is stashed
  std::string solution_file = libMesh_inputfile( "ExactSolution/solution_file", "DIE!" );
  es_ref.read( solution_file );

  exact_sol.attach_reference_solution( &es_ref );
  
  // Compute error and get it in various norms
  exact_sol.compute_error("GRINS", "u");
  exact_sol.compute_error("GRINS", "v");

  if( (es->get_mesh()).mesh_dimension() == 3 )
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

  if( (es->get_mesh()).mesh_dimension() == 3 )
    {
      w_l2error = exact_sol.l2_error("GRINS", "w");
      w_h1error = exact_sol.h1_error("GRINS", "w");
    }

  int return_flag = 0;

  // This is the tolerance of the iterative linear solver so
  // it's unreasonable to expect anything better than this.
  double tol = 5.0e-13;
  
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
