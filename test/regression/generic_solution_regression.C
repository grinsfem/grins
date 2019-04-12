//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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
#include "grins/runner.h"
#include "grins/mesh_builder.h"
#include "grins/multiphysics_sys.h"

//libMesh
#include "libmesh/exact_solution.h"

void test_error_norm( libMesh::ExactSolution& exact_sol,
                      const std::string& system_name,
                      const std::string& var,
                      const std::string& norm,
                      const double tol,
                      int& return_flag );

int main(int argc, char* argv[])
{
  GRINS::Runner grins(argc,argv);

  const GetPot & command_line = grins.get_command_line();

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

  const GetPot & inputfile = grins.get_input_file();

  // Don't flag our command-line-specific variables as UFOs later
  inputfile.have_variable("soln-data");
  inputfile.have_variable("vars");
  inputfile.have_variable("norms");
  inputfile.have_variable("tol");
  inputfile.have_variable("qois");

  // Initialize Simulation
  grins.init();

  // Do solve here
  grins.run();

  GRINS::Simulation & sim = grins.get_simulation();

  // Get equation systems to create ExactSolution object
  std::shared_ptr<libMesh::EquationSystems> es = sim.get_equation_system();

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

  const std::string& system_name = sim.get_multiphysics_system_name();

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

  // Now grab any "gold" QoI values
  unsigned int n_qoi_vals = command_line.vector_variable_size("qois");
  for( unsigned int n = 0; n < n_qoi_vals; n++ )
    {

      std::cout << "==========================================================" << std::endl
                << "Checking qoi " << n << " with tol " << tol << "...";

      libMesh::Number gold_qoi = command_line("qois", libMesh::Number(0), n);

      libMesh::Number computed_qoi =
        sim.get_multiphysics_system()->qoi[n];

      double error = computed_qoi - gold_qoi;

      if (std::abs(error) > tol)
        {
          std::cerr << "Tolerance exceeded for generic regression test!" << std::endl
                    << "tolerance     = " << tol << std::endl
                    << "error         = " << error << std::endl
                    << "qoi index     = " << n << std::endl;
          return_flag = 1;
        }
      else
        std::cout << "PASSED!" << std::endl
                  << "==========================================================" << std::endl;

    }


  return return_flag;
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
