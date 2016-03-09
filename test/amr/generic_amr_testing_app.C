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

  if( !command_line.have_variable("test_data_prefix") )
    {
      std::cerr << "ERROR: Must specify solution test data filename prefix on command line with test_data_prefix=<string>."
                << std::endl;
      exit(1);
    }

  if( !command_line.have_variable("mesh_data_prefix") )
    {
      std::cerr << "ERROR: Must specify solution mesh data filename prefix on command line with test_data_prefix=<string>."
                << std::endl;
      exit(1);
    }

  if( !command_line.have_variable("gold_data_prefix") )
    {
      std::cerr << "ERROR: Must specify solution gold data filename prefix on command line with gold_data_prefix=<string>."
                << std::endl;
      exit(1);
    }

  if( !command_line.have_variable("gold_mesh_prefix") )
    {
      std::cerr << "ERROR: Must specify gold mesh data filename prefix on command line with gold_mesh_prefix=<string>."
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
      std::cerr << "ERROR: Must specify test tolerance on command line with tol=<double>"
                << std::endl;
      exit(1);
    }

  if( !command_line.have_variable("n_steps") )
    {
      std::cerr << "ERROR: Must specify n_steps on command line with n_steps=<int>"
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

  std::string test_data_prefix = command_line("test_data_prefix", "DIE!" );
  std::string mesh_data_prefix = command_line("mesh_data_prefix", "DIE!" );

  std::string gold_data_prefix = command_line("gold_data_prefix", "DIE!" );
  std::string gold_mesh_prefix = command_line("gold_mesh_prefix", "DIE!" );

  unsigned int n_steps = command_line("n_steps", 0);

  int return_flag = 0;

  std::string system_name = "GRINS-TEST";

  for( unsigned int s = 0; s < n_steps; s++ )
    {
      std::stringstream step_string;
      step_string << s;

      //FIXME: Need to support different input formats for meshes
      std::string mesh_filename = mesh_data_prefix+"."+step_string.str()+"_mesh.xda";
      {
        std::ifstream i( mesh_filename.c_str());
        if (!i)
          {
            std::cerr << "Error: Could not read from mesh_data file "
                      << mesh_filename << std::endl;
            exit(1);
          }
      }

      libMesh::Mesh mesh(libmesh_init.comm());
      mesh.read(mesh_filename);

      // This needs to match the read counter-part in GRINS::Simulation
      //FIXME: Need to support different input formats for restarts
      std::string test_data = test_data_prefix+"."+step_string.str()+".xdr";

      {
        std::ifstream i(test_data.c_str());
        if (!i)
          {
            std::cerr << "Error: Could not read from test_data file "
                      << test_data << std::endl;
            exit(1);
          }
      }

      libMesh::EquationSystems es(mesh);
      es.read(test_data,
              GRINSEnums::DECODE,
              libMesh::EquationSystems::READ_HEADER |
              libMesh::EquationSystems::READ_DATA |
              libMesh::EquationSystems::READ_ADDITIONAL_DATA);

      // Now grab the "gold" mesh and data
      //FIXME: Need to support different input formats for meshes
      std::string gold_mesh_filename = gold_mesh_prefix+"_mesh."+step_string.str()+".exo";
      {
        std::ifstream i( gold_mesh_filename.c_str());
        if (!i)
          {
            std::cerr << "Error: Could not read from gold_mesh file "
                      << gold_mesh_filename << std::endl;
            exit(1);
          }
      }

      libMesh::Mesh gold_mesh(libmesh_init.comm());
      gold_mesh.read(gold_mesh_filename);

      // This needs to match the read counter-part in GRINS::Simulation
      //FIXME: Need to support different input formats for restarts
      std::string gold_data = gold_data_prefix+"."+step_string.str()+".xdr";

      {
        std::ifstream i(gold_data.c_str());
        if (!i)
          {
            std::cerr << "Error: Could not read from test_data file "
                      << gold_data << std::endl;
            exit(1);
          }
      }

      libMesh::EquationSystems es_gold(gold_mesh);
      es_gold.read(gold_data,
                   GRINSEnums::DECODE,
                   libMesh::EquationSystems::READ_HEADER |
                   libMesh::EquationSystems::READ_DATA |
                   libMesh::EquationSystems::READ_ADDITIONAL_DATA);

      // Create Exact solution object and attach "gold" EquationSystems
      libMesh::ExactSolution exact_sol(es);

      exact_sol.attach_reference_solution( &es_gold );

      // Now compute error for each variable
      for( unsigned int v = 0; v < n_vars; v++ )
        {
          exact_sol.compute_error(system_name, vars[v]);
        }

      double tol = command_line("tol", 1.0e-10);

      return_flag = 0;
      int test_flag = 0;

      // We're just diffing against gold data, so error should be 0
      double exact_error = 0.0;

      // Now test error for each variable, for each norm
      for( unsigned int v = 0; v < n_vars; v++ )
        {
          for( unsigned int n = 0; n < n_norms; n++ )
            {
              test_flag = test_error_norm( exact_sol, system_name, vars[v], norms[n], exact_error, tol );

              if( test_flag == 1 )
                return_flag = 1;
            }
        }

      // Break out if one of the tests failed
      if( return_flag == 1 )
        break;
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
      std::cerr << "ERROR: Invalid norm " << norm << std::endl;
      exit(1);
    }

  if( std::abs(error-exact_error) > tol )
    {
      return_flag = 1;

      std::cerr << "Tolerance exceeded for generic regression test!" << std::endl
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
