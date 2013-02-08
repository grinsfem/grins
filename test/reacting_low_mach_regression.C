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
#include "grins_config.h"

#include <iostream>

// GRINS
#include "grins/simulation.h"
#include "grins/simulation_builder.h"

// GRVY
#ifdef HAVE_GRVY
#include "grvy.h"
#endif

// libMesh
#include "libmesh/parallel.h"
#include "libmesh/exact_solution.h"

// Function for getting initial temperature field
Real initial_values( const Point& p, const Parameters &params, 
		     const std::string& system_name, const std::string& unknown_name );

int run( int argc, char* argv[], const GetPot& input );

int main(int argc, char* argv[])
{
#ifdef USE_GRVY_TIMERS
  GRVY::GRVY_Timer_Class grvy_timer;
  grvy_timer.Init("GRINS Timer");
#endif

  // Check command line count.
  if( argc < 3 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify libMesh input file and exact solution file." << std::endl;
      exit(1); // TODO: something more sophisticated for parallel runs?
    }

  // libMesh input file should be first argument
  std::string libMesh_input_filename = argv[1];
  
  // Create our GetPot object.
  GetPot libMesh_inputfile( libMesh_input_filename );

  int return_flag = 0;

  // If we have Cantera, always run
#ifdef GRINS_HAVE_CANTERA
  return_flag = run(argc,argv,libMesh_inputfile);
#else
  if( std::string( libMesh_inputfile("Physics/ReactingLowMachNavierStokes/chemistry_library", "DIE!") ) == std::string("cantera") )
    {
      return_flag = 77;
    }
  else
    {
      return_flag = run(argc,argv,libMesh_inputfile);
    }
#endif

  return return_flag;
}

int run( int argc, char* argv[], const GetPot& input )
{
#ifdef USE_GRVY_TIMERS
  grvy_timer.BeginTimer("Initialize Solver");
#endif

  // Initialize libMesh library.
  LibMeshInit libmesh_init(argc, argv);
 
  GRINS::SimulationBuilder sim_builder;

  GRINS::Simulation grins( input,
			   sim_builder );

  //FIXME: We need to move this to within the Simulation object somehow...
  std::string restart_file = input( "restart-options/restart_file", "none" );

  if( restart_file == "none" )
    {
      // Asssign initial temperature value
      std::string system_name = input( "screen-options/system_name", "GRINS" );
      std::tr1::shared_ptr<libMesh::EquationSystems> es = grins.get_equation_system();
      const libMesh::System& system = es->get_system(system_name);
      
      Parameters &params = es->parameters;

      Real& w_N2 = params.set<Real>( "w_N2" );
      w_N2 = input( "Physics/ReactingLowMachNavierStokes/bound_species_3", 0.0, 0 );
      
      Real& w_N = params.set<Real>( "w_N" );
      w_N = input( "Physics/ReactingLowMachNavierStokes/bound_species_3", 0.0, 1 );

      system.project_solution( initial_values, NULL, params );
    }

#ifdef USE_GRVY_TIMERS
  grvy_timer.EndTimer("Initialize Solver");

  // Attach GRVY timer to solver
  grins.attach_grvy_timer( &grvy_timer );
#endif

  grins.run();

#ifdef USE_GRVY_TIMERS
  grvy_timer.Finalize();
 
  if( Parallel::Communicator_World.rank() == 0 ) grvy_timer.Summarize();
#endif

  // Get equation systems to create ExactSolution object
  std::tr1::shared_ptr<EquationSystems> es = grins.get_equation_system();

  //es->write("foobar.xdr");

  // Create Exact solution object and attach exact solution quantities
  ExactSolution exact_sol(*es);
  
  EquationSystems es_ref( es->get_mesh() );

  // Filename of file where comparison solution is stashed
  std::string solution_file = std::string(argv[2]);
  es_ref.read( solution_file );

  exact_sol.attach_reference_solution( &es_ref );
  
  // Compute error and get it in various norms
  exact_sol.compute_error("GRINS", "u");
  exact_sol.compute_error("GRINS", "v");

  if( (es->get_mesh()).mesh_dimension() == 3 )
    exact_sol.compute_error("GRINS", "w");

  exact_sol.compute_error("GRINS", "p");
  exact_sol.compute_error("GRINS", "T");
  exact_sol.compute_error("GRINS", "w_N2");
  exact_sol.compute_error("GRINS", "w_N");

  double u_l2error = exact_sol.l2_error("GRINS", "u");
  double u_h1error = exact_sol.h1_error("GRINS", "u");

  double v_l2error = exact_sol.l2_error("GRINS", "v");
  double v_h1error = exact_sol.h1_error("GRINS", "v");

  double p_l2error = exact_sol.l2_error("GRINS", "p");
  double p_h1error = exact_sol.h1_error("GRINS", "p");

  double T_l2error = exact_sol.l2_error("GRINS", "T");
  double T_h1error = exact_sol.h1_error("GRINS", "T");

  double wN_l2error = exact_sol.l2_error("GRINS", "w_N");
  double wN_h1error = exact_sol.h1_error("GRINS", "w_N");

  double wN2_l2error = exact_sol.l2_error("GRINS", "w_N2");
  double wN2_h1error = exact_sol.h1_error("GRINS", "w_N2");
  
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
  double tol = 1.0e-10;
  
  if( u_l2error > tol   || u_h1error > tol   ||
      v_l2error > tol   || v_h1error > tol   ||
      w_l2error > tol   || w_h1error > tol   ||
      p_l2error > tol   || p_h1error > tol   ||
      T_l2error > tol   || T_h1error > tol   ||
      wN_l2error > tol  || wN_h1error > tol  ||
      wN2_l2error > tol || wN2_h1error > tol    )
    {
      return_flag = 1;

      std::cout << "Tolerance exceeded for thermally driven flow test." << std::endl
		<< "tolerance     = " << tol << std::endl
		<< "u l2 error    = " << u_l2error << std::endl
		<< "u h1 error    = " << u_h1error << std::endl
		<< "v l2 error    = " << v_l2error << std::endl
		<< "v h1 error    = " << v_h1error << std::endl
		<< "w l2 error    = " << w_l2error << std::endl
		<< "w h1 error    = " << w_h1error << std::endl
		<< "p l2 error    = " << p_l2error << std::endl
		<< "p h1 error    = " << p_h1error << std::endl
		<< "T l2 error    = " << T_l2error << std::endl
		<< "T h1 error    = " << T_h1error << std::endl
		<< "w_N l2 error  = " << wN_l2error << std::endl
		<< "w_N h1 error  = " << wN_h1error << std::endl
		<< "w_N2 l2 error = " << wN2_l2error << std::endl
		<< "w_N2 h1 error = " << wN2_h1error << std::endl;
    }

  return return_flag;
}

Real initial_values( const Point& p, const Parameters &params, 
		     const std::string& , const std::string& unknown_name )
{
  Real value = 0.0;

  if( unknown_name == "w_N2" )
    value = params.get<Real>("w_N2");

  else if( unknown_name == "w_N" )
    value = params.get<Real>("w_N");

  else if( unknown_name == "T" )
    value = 300;

  else if( unknown_name == "u" )
    {
      const Real y = p(1);
      value = 1.0*(-y*y+1);
    }

  else
    value = 0.0;

  return value;
}
