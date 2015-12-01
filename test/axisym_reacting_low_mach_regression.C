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
#include "grins/simulation.h"
#include "grins/simulation_builder.h"

// libMesh
#include "libmesh/parallel.h"
#include "libmesh/exact_solution.h"

// Function for getting initial temperature field
libMesh::Real
initial_values( const libMesh::Point& p, const libMesh::Parameters &params,
		const std::string& system_name, const std::string& unknown_name );

int run( int argc, char* argv[], const GetPot& input, GetPot& command_line );

void test_error_norm( libMesh::ExactSolution& exact_sol,
                      const std::string& system_name,
                      const std::string& var,
                      const std::string& norm,
                      const double tol,
                      int& return_flag );

int main(int argc, char* argv[])
{
  GetPot command_line(argc,argv);

  if( !command_line.have_variable("input") )
    {
      std::cerr << "ERROR: Must specify input file on command line with input=<file>." << std::endl;
      exit(1);
    }

  if( !command_line.have_variable("soln-data") )
    {
      std::cerr << "ERROR: Must specify solution data on command line with soln-data=<file>." << std::endl;
      exit(1);
    }

  if( !command_line.have_variable("vars") )
    {
      std::cerr << "ERROR: Must specify variables on command line with vars='var1 var2'" << std::endl;
      exit(1);
    }

  if( !command_line.have_variable("norms") )
    {
      std::cerr << "ERROR: Must specify error norms on command line with norms='L2 H1'" << std::endl;
      exit(1);
    }

  if( !command_line.have_variable("tol") )
    {
      std::cerr << "ERROR: Must specify test tolerance on command line with tol=<tol>" << std::endl;
      exit(1);
    }

  // libMesh input file should be first argument
  std::string libMesh_input_filename = command_line("input", "DIE!");

  {
    std::ifstream i(libMesh_input_filename.c_str());
    if (!i)
      {
        std::cerr << "Error: Could not read from libMesh input file "
                << libMesh_input_filename << std::endl;
        exit(1);
      }
  }

  // Create our GetPot object.
  GetPot libMesh_inputfile( libMesh_input_filename );

  int return_flag = 0;

  std::string chem_lib = libMesh_inputfile("Physics/ReactingLowMachNavierStokes/thermochemistry_library", "DIE!");

  if( chem_lib == std::string("cantera") )
    {
#ifdef GRINS_HAVE_CANTERA
      return_flag = run(argc,argv,libMesh_inputfile,command_line);
#else
      return_flag = 77;
#endif
    }
  else if( chem_lib == std::string("antioch") )
    {
#ifdef GRINS_HAVE_ANTIOCH
      return_flag = run(argc,argv,libMesh_inputfile,command_line);
#else
      return_flag = 77;
#endif
    }
  else
    {
      return_flag = 1;
    }

  return return_flag;
}

int run( int argc, char* argv[], const GetPot& input, GetPot& command_line )
{
  // Initialize libMesh library.
  libMesh::LibMeshInit libmesh_init(argc, argv);

  GRINS::SimulationBuilder sim_builder;

  GRINS::Simulation grins( input,
			   sim_builder,
                           libmesh_init.comm() );

  //FIXME: We need to move this to within the Simulation object somehow...
  std::string restart_file = input( "restart-options/restart_file", "none" );

  if( restart_file == "none" )
    {
      // Asssign initial temperature value
      std::string system_name = input( "screen-options/system_name", "GRINS" );
      GRINS::SharedPtr<libMesh::EquationSystems> es = grins.get_equation_system();
      const libMesh::System& system = es->get_system(system_name);

      libMesh::Parameters &params = es->parameters;

      libMesh::Real& w_N2 = params.set<libMesh::Real>( "w_N2" );
      w_N2 = input( "Physics/ReactingLowMachNavierStokes/bound_species_0", 0.0, 0.0 );

      libMesh::Real& w_N = params.set<libMesh::Real>( "w_N" );
      w_N = input( "Physics/ReactingLowMachNavierStokes/bound_species_0", 0.0, 1.0 );

      system.project_solution( initial_values, NULL, params );
    }

  grins.run();

  // Get equation systems to create ExactSolution object
  GRINS::SharedPtr<libMesh::EquationSystems> es = grins.get_equation_system();

   // Create Exact solution object and attach exact solution quantities
  libMesh::ExactSolution exact_sol(*es);

  libMesh::EquationSystems es_ref( es->get_mesh() );

  // Filename of file where comparison solution is stashed
  std::string solution_file = command_line("soln-data", "DIE!");
  es_ref.read( solution_file );

  exact_sol.attach_reference_solution( &es_ref );

  // Now grab the variables for which we want to compare
  unsigned int n_vars = command_line.vector_variable_size("vars");
  std::vector<std::string> vars(n_vars);
  for( unsigned int v = 0; v < n_vars; v++ )
    {
      vars[v] = command_line("vars", "DIE!", v);
    }

  // Now grab the norms to compute for each variable error
  unsigned int n_norms = command_line.vector_variable_size("norms");
  std::vector<std::string> norms(n_norms);
  for( unsigned int n = 0; n < n_norms; n++ )
    {
      norms[n] = command_line("norms", "DIE!", n);
      if( norms[n] != std::string("L2") &&
          norms[n] != std::string("H1") )
        {
          std::cerr << "ERROR: Invalid norm input " << norms[n] << std::endl
                    << "       Valid values are: L2" << std::endl
                    << "                         H1" << std::endl;
        }
    }

  const std::string& system_name = grins.get_multiphysics_system_name();

  // Now compute error for each variable
  for( unsigned int v = 0; v < n_vars; v++ )
    {
      exact_sol.compute_error(system_name, vars[v]);
    }

  int return_flag = 0;

  double tol = command_line("tol", 1.0e-10);

  // Now test error for each variable, for each norm
  for( unsigned int v = 0; v < n_vars; v++ )
    {
      for( unsigned int n = 0; n < n_norms; n++ )
        {
          test_error_norm( exact_sol, system_name, vars[v], norms[n], tol, return_flag );
        }
    }

  return return_flag;
}

libMesh::Real
initial_values( const libMesh::Point& p, const libMesh::Parameters &params,
		const std::string& , const std::string& unknown_name )
{
  libMesh::Real value = 0.0;

  if( unknown_name == "w_N2" )
    value = params.get<libMesh::Real>("w_N2");

  else if( unknown_name == "w_N" )
    value = params.get<libMesh::Real>("w_N");

  else if( unknown_name == "T" )
    value = 300;

  else if( unknown_name == "u" )
    {
      const libMesh::Real x = p(0);
      value = 1.0*(-x*x+1);
    }

  else
    value = 0.0;

  return value;
}

void test_error_norm( libMesh::ExactSolution& exact_sol,
                      const std::string& system_name,
                      const std::string& var,
                      const std::string& norm,
                      const double tol,
                      int& return_flag )
{
  // We don't set return_flag unless we are setting it 1
  // since this function gets called multiple times and we don't
  // want to overwrite a previous "fail" (return_flag = 1) with
  // a "pass" (return_flag = 0)

  double error = 0.0;

  std::cout << "==========================================================" << std::endl
            << "Checking variable " << var << " using error norm " << norm << " with tol " << tol << "...";

  if( norm == std::string("L2") )
    {
      error = exact_sol.l2_error(system_name, var);
    }
  else if( norm == std::string("H1") )
    {
      error = exact_sol.h1_error(system_name, var);
    }
  else
    {
      std::cerr << "ERROR: Invalid norm " << norm << std::endl;
      exit(1);
    }

  if( error > tol )
    {
      return_flag = 1;

      std::cerr << "Tolerance exceeded for generic regression test!" << std::endl
                << "tolerance     = " << tol << std::endl
                << "norm of error = " << error << std::endl
                << "norm type     = " << norm << std::endl
                << "var           = " << var << std::endl;
    }
  else
    {
      std::cout << "PASSED!" << std::endl
                << "==========================================================" << std::endl;
    }

  return;
}
