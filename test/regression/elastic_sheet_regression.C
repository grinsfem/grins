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

#include <iostream>

// GRINS
#include "grins/runner.h"

// libMesh
#include "libmesh/parallel.h"
#include "libmesh/exact_solution.h"

// Function for getting initial temperature field
libMesh::Real
initial_values( const libMesh::Point& p, const libMesh::Parameters &params,
                const std::string& system_name, const std::string& unknown_name );

int run( char* argv[], GRINS::Runner & grins );

int main(int argc, char* argv[])
{
  // Check command line count.
  if( argc < 3 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify libMesh input file and exact solution file." << std::endl;
      exit(1); // TODO: something more sophisticated for parallel runs?
    }

  GRINS::Runner grins(argc,argv);

  int return_flag = 0;

  return_flag = run( argv, grins );

  return return_flag;
}

int run( char* argv[], GRINS::Runner & grins )
{
  grins.init();
  grins.run();

  GRINS::Simulation & sim = grins.get_simulation();

  // Get equation systems to create ExactSolution object
  std::shared_ptr<libMesh::EquationSystems>
    es = sim.get_equation_system();

  //es->write("foobar.xdr");

  // Create Exact solution object and attach exact solution quantities
  libMesh::ExactSolution exact_sol(*es);

  libMesh::EquationSystems es_ref( es->get_mesh() );

  // Filename of file where comparison solution is stashed
  std::string solution_file = std::string(argv[2]);
  es_ref.read( solution_file );

  exact_sol.attach_reference_solution( &es_ref );

  // Compute error and get it in various norms
  exact_sol.compute_error("StretchedElasticSheet", "u");
  exact_sol.compute_error("StretchedElasticSheet", "v");
  exact_sol.compute_error("StretchedElasticSheet", "w");

  double u_l2error = exact_sol.l2_error("StretchedElasticSheet", "u");
  double u_h1error = exact_sol.h1_error("StretchedElasticSheet", "u");

  double v_l2error = exact_sol.l2_error("StretchedElasticSheet", "v");
  double v_h1error = exact_sol.h1_error("StretchedElasticSheet", "v");

  double w_l2error = exact_sol.l2_error("StretchedElasticSheet", "w");
  double w_h1error = exact_sol.h1_error("StretchedElasticSheet", "w");

  int return_flag = 0;

  double tol = 5.0e-8;

  if( u_l2error > tol   || u_h1error > tol   ||
      v_l2error > tol   || v_h1error > tol   ||
      w_l2error > tol   || w_h1error > tol     )
    {
      return_flag = 1;

      std::cout << "Tolerance exceeded for thermally driven flow test." << std::endl
                << "tolerance     = " << tol << std::endl
                << "u l2 error    = " << u_l2error << std::endl
                << "u h1 error    = " << u_h1error << std::endl
                << "v l2 error    = " << v_l2error << std::endl
                << "v h1 error    = " << v_h1error << std::endl
                << "w l2 error    = " << w_l2error << std::endl
                << "w h1 error    = " << w_h1error << std::endl;
    }

  return return_flag;
}
