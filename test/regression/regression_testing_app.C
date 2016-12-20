//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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
#include "grins/mesh_builder.h"
#include "grins/simulation.h"
#include "grins/simulation_builder.h"
#include "grins/multiphysics_sys.h"
#include "grins/string_utils.h"

//libMesh
#include "libmesh/exact_solution.h"

void test_error_norm( libMesh::ExactSolution& exact_sol,
                      const std::string& system_name,
                      const std::string& var,
                      const std::string& norm,
                      const double tol,
                      int& return_flag );

void parse_qoi_data( const std::vector<std::string> & qoi_names,
                     const GetPot & command_line,
                     std::vector<libMesh::Real> & qoi_data_values,
                     std::vector<libMesh::Real> & qoi_gold_values );

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

  if( !command_line.have_variable("gold-data") )
    {
      std::cerr << "ERROR: Must specify gold data on command line with gold-data=<file>." << std::endl;
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

  if( command_line.have_variable("gold-qoi-names") && !command_line.have_variable("qoi-data") )
    {
      std::cerr << "ERROR: Must specify QoI data file on command line with qoi-data=<file>" << std::endl;
      exit(1);
    }

  if( command_line.have_variable("gold-qoi-names") && !command_line.have_variable("gold-qoi-values") )
    {
      std::cerr << std::string("ERROR: Must specify QoI gold values on command line with gold-qoi-values='qoi1 qoi2'")+
        std::string(" in the same order as given in gold-qoi-names") << std::endl;
      exit(1);
    }

  if( command_line.have_variable("gold-qoi-names") )
    {
      // Make sure there are the same number of values as names
      unsigned int n_qoi_names = command_line.vector_variable_size("gold-qoi-names");
      unsigned int n_qoi_values = command_line.vector_variable_size("gold-qoi-values");

      if( n_qoi_names != n_qoi_values )
        {
          std::cerr << "ERROR: Must have the same number of gold-qoi-names and gold-qoi-values!" << std::endl;
          exit(1);
        }


    }

  // Now grab the variables for which we want to compare
  unsigned int n_vars = command_line.vector_variable_size("vars");
  std::vector<std::string> vars(n_vars);
  for( unsigned int v = 0; v < n_vars; v++ )
    vars[v] = command_line("vars", "DIE!", v);

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

  // Parse QoI information
  std::vector<std::string> qoi_names;
  std::vector<libMesh::Real> qoi_data_values;
  std::vector<libMesh::Real> qoi_gold_values;

  if( command_line.have_variable("gold-qoi-names") )
    {
      unsigned int n_qoi_names = command_line.vector_variable_size("gold-qoi-names");
      qoi_names.reserve(n_qoi_names);
      for( unsigned int i = 0; i < n_qoi_names; i++ )
        qoi_names.push_back( command_line("gold-qoi-names", "DIE!") );


      parse_qoi_data( qoi_names, command_line, qoi_data_values, qoi_gold_values );
    }

  std::string input_filename = command_line("input", "DIE!");

  {
    std::ifstream i(input_filename.c_str());
    if (!i)
      {
        std::cerr << "Error: Could not read from libMesh input file "
                << input_filename << std::endl;
        exit(1);
      }
  }

  // Create our GetPot object.
  GetPot input( input_filename );

  // But allow command line options to override the file
  input.parse_command_line(argc, argv);

  // Don't flag our command-line-specific variables as UFOs later
  input.have_variable("input");
  input.have_variable("soln-data");
  input.have_variable("gold-data");
  input.have_variable("vars");
  input.have_variable("norms");
  input.have_variable("tol");
  input.have_variable("gold-qoi-names");
  input.have_variable("qoi-data");
  input.have_variable("gold-qoi-values");

  // Initialize libMesh library.
  libMesh::LibMeshInit libmesh_init(argc, argv);

  GRINS::MeshBuilder mesh_builder;
  GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh = mesh_builder.build(input,libmesh_init.comm());

  libMesh::EquationSystems es(*mesh);

  // This needs to match the read counter-part in GRINS::Simulation
  std::string soln_data = command_line("soln-data", "DIE!" );

  {
    std::ifstream i(soln_data.c_str());
    if (!i)
      {
        std::cerr << "Error: Could not read from soln_data file "
                  << soln_data << std::endl;
        exit(1);
      }
  }

  es.read(soln_data,
          libMesh::EquationSystems::READ_HEADER |
          libMesh::EquationSystems::READ_DATA |
          libMesh::EquationSystems::READ_ADDITIONAL_DATA);

  std::string system_name = "GRINS-TEST";

  // Now build up the gold data
  libMesh::EquationSystems es_ref(*mesh);

  // Filename of file where comparison solution is stashed
  std::string gold_file = command_line("gold-data", "DIE!");
  es_ref.read( gold_file );

  // Create Exact solution object and attach exact solution quantities
  libMesh::ExactSolution exact_sol(es);
  exact_sol.attach_reference_solution( &es_ref );

  // Now compute error for each variable
  for( unsigned int v = 0; v < n_vars; v++ )
    exact_sol.compute_error(system_name, vars[v]);

  int return_flag = 0;

  double tol = command_line("tol", 1.0e-10);

  // Now test error for each variable, for each norm
  for( unsigned int v = 0; v < n_vars; v++ )
    for( unsigned int n = 0; n < n_norms; n++ )
      test_error_norm( exact_sol, system_name, vars[v], norms[n], tol, return_flag );

  // Compare QoI values
  unsigned int n_qois = qoi_names.size();

  for( unsigned int n = 0; n < n_qois; n++ )
    {

      std::cout << "==========================================================" << std::endl
                << "Checking qoi " << qoi_names[n] << " with tol " << tol << "...";

      libMesh::Number gold_qoi = qoi_gold_values[n];

      libMesh::Number computed_qoi = qoi_data_values[n];

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
}

void parse_qoi_data( const std::vector<std::string> & qoi_names,
                     const GetPot & command_line,
                     std::vector<libMesh::Real> & qoi_data_values,
                     std::vector<libMesh::Real> & qoi_gold_values )
{
  std::string qoi_data_filename = command_line("qoi-data", "nofile!");

  std::ifstream qoi_input( qoi_data_filename );
  if( !qoi_input.good() )
    {
      std::cerr << "ERROR: Could not open "+qoi_data_filename << std::endl;
      exit(1);
    }

  unsigned int n_qois = qoi_names.size();
  qoi_data_values.resize(n_qois);
  qoi_gold_values.resize(n_qois);

  // Grab the first line of the file
  std::string line;
  std::getline( qoi_input, line );

  while( !qoi_input.eof() )
    {
      std::vector<std::string> tokens;
      GRINS::StringUtilities::split_string( line, " =", tokens );

      std::string qoi_name = tokens[0];
      libMesh::Real qoi_value = GRINS::StringUtilities::string_to_T<libMesh::Real>(tokens[1]);

      // Figure out the index in the qoi_names vector
      int idx = std::find(qoi_names.begin(), qoi_names.end(), qoi_name) - qoi_names.begin();
      if( idx < (int)qoi_names.size() )
        qoi_data_values[idx] = qoi_value;

      // Grab the next line. If we hit the end of the file, this will flip the eofbit.
      std::getline( qoi_input, line );
    }

  qoi_input.close();
}
