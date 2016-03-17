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

// GRINS
#include "grins/grins_enums.h"
#include "grins/mesh_builder.h"
#include "grins/simulation.h"
#include "grins/simulation_builder.h"
#include "grins/multiphysics_sys.h"

//libMesh
#include "libmesh/exact_solution.h"
#include "libmesh/parsed_function.h"
#include "libmesh/composite_function.h"

int test_error_norm( libMesh::ExactSolution& exact_sol,
                     const std::string& system_name,
                     const std::string& var,
                     const std::string& norm,
                     const libMesh::Real exact_error,
                     const double tol );

int main(int argc, char* argv[])
{
  GetPot command_line(argc,argv);

  if( !command_line.have_variable("input") )
    {
      std::cerr << "ERROR: Must specify input file on command line with input=<file>."
                << std::endl;
      exit(1);
    }

  if( !command_line.have_variable("test_data") )
    {
      std::cerr << "ERROR: Must specify solution test data on command line with test_data=<file>."
                << std::endl;
      exit(1);
    }

  if( !command_line.have_variable("vars") )
    {
      std::cerr << "ERROR: Must specify variables on command line with vars='var1 var2'"
                << std::endl;
      exit(1);
    }

  if( !command_line.have_variable("norms") )
    {
      std::cerr << "ERROR: Must specify error norms on command line with norms='L2 H1'"
                << std::endl;
      exit(1);
    }

  if( !command_line.have_variable("tol") )
    {
      std::cerr << "ERROR: Must specify test tolerance on command line with tol=<tol>"
                << std::endl;
      exit(1);
    }



  // Grab the variables for which we want to compare the errors
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
      if( norms[n] != std::string("L2") )
        {
          std::cerr << "ERROR: Invalid norm input " << norms[n] << std::endl
                    << "       Valid values are: L2" << std::endl;
        }
    }

  // Now that we have the norms, we need to grab the error values for each var, for each norm from the CLI
  std::vector<std::vector<libMesh::Real> > error_values;
  error_values.resize(n_vars);

  for( unsigned int v = 0; v < n_vars; v++ )
    {
      std::string var = vars[v];

      error_values[v].resize(n_norms);

      for( unsigned int n = 0; n < n_norms; n++ )
        {
          std::string norm_cli = var+"_"+norms[n]+"_error";

          if( !command_line.have_variable(norm_cli) )
            {
              std::cerr << "ERROR: Must specify "+norms[n]+" errors on command line with "+
                var+"_"+norms[n]+"_error='<error>'" << std::endl;
              exit(1);
            }

          error_values[v][n] = command_line( norm_cli, -1.0 );

          if( error_values[v][n] < 0.0 )
            {
              std::cerr << "ERROR: Invalid value of error --- norms should be positive!" << std::endl
                        << "       Found error value: " << error_values[v][n] << std::endl;
              exit(1);
            }
        }
    }

  // Grab the input filename
  std::string input_filename = command_line("input", "DIE!");

  // Check that the input file is actually there
  {
    std::ifstream i(input_filename.c_str());
    if (!i)
      {
        std::cerr << "Error: Could not read from input file "
                  << input_filename << std::endl;
        exit(1);
      }
  }

  // Create our GetPot object.
  GetPot input( input_filename );

  // Initialize libMesh library
  libMesh::LibMeshInit libmesh_init(argc, argv);

  GRINS::MeshBuilder mesh_builder;
  GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh = mesh_builder.build(input,libmesh_init.comm());

  libMesh::EquationSystems es(*mesh);

  // This needs to match the read counter-part in GRINS::Simulation
  std::string test_data = command_line("test_data", "DIE!" );
  
  {
    std::ifstream i(test_data.c_str());
    if (!i)
      {
        std::cerr << "Error: Could not read from test_data file "
                  << test_data << std::endl;
        exit(1);
      }
  }

  es.read(test_data,
          libMesh::EquationSystems::READ_HEADER |
          libMesh::EquationSystems::READ_DATA |
          libMesh::EquationSystems::READ_ADDITIONAL_DATA);

  std::string system_name = "GRINS-TEST";
  libMesh::System& system = es.get_system(system_name);

  // Create Exact solution object and attach exact solution quantities
  libMesh::ExactSolution exact_sol(es);

  // Grab the exact solution expression for each variable
  // We'll directly construct the ParsedFunction for each of the variable exact solutions provided
  // They are all packed together in the CompositeFunction
  // We need to do this here because we need the variable index from the equation system
  {
    libMesh::CompositeFunction<libMesh::Real> exact_sols;

    for( unsigned int v = 0; v < n_vars; v++ )
      {
        std::string var = vars[v];
        std::string soln_cli = var+"_exact_soln";

        if( !command_line.have_variable(soln_cli) )
          {
            std::cerr << "ERROR: Must specify the exact solution for the variable"+var+"on" << std::endl
                      << "       the command line with "+soln_cli+"=<expression>"
                      << std::endl;
            exit(1);
          }

        unsigned int var_index = system.variable_number(var);

        std::string parsed_soln = command_line(soln_cli,"NULL");

        std::vector<unsigned int> index_map(1,var_index);

        exact_sols.attach_subfunction( libMesh::ParsedFunction<libMesh::Real>(parsed_soln), index_map );
      }

    // This is assuming there's only 1 system
    exact_sol.attach_exact_value(0, &exact_sols);
  }

  // Now compute error for each variable
  for( unsigned int v = 0; v < n_vars; v++ )
    {
      exact_sol.compute_error(system_name, vars[v]);
    }

  double tol = command_line("tol", 1.0e-10);

  int return_flag = 0;
  int test_flag = 0;

  // Now test error for each variable, for each norm
  for( unsigned int v = 0; v < n_vars; v++ )
    {
      for( unsigned int n = 0; n < n_norms; n++ )
        {
          test_flag = test_error_norm( exact_sol, system_name, vars[v], norms[n], error_values[v][n], tol );

          if( test_flag == 1 )
            return_flag = 1;
        }
    }

  return return_flag;
}

int test_error_norm( libMesh::ExactSolution& exact_sol,
                      const std::string& system_name,
                      const std::string& var,
                      const std::string& norm,
                      const libMesh::Real exact_error,
                      const double tol )
{
  int return_flag = 0;

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
      std::cout << "ERROR: Invalid norm " << norm << std::endl;
      exit(1);
    }

  if( std::abs(error-exact_error) > tol )
    {
      return_flag = 1;

      std::cout << "Tolerance exceeded for generic regression test!" << std::endl
                << "tolerance     = " << tol << std::endl
                << "norm of error = " << error << std::endl
                << "exact error   = " << exact_error << std::endl
                << "norm type     = " << norm << std::endl
                << "var           = " << var << std::endl;
    }
  else
    {
      std::cout << "PASSED!" << std::endl
                << "==========================================================" << std::endl;
    }

  return return_flag;
}
