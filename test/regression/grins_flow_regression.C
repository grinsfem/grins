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


// C++
#include <iostream>

// GRINS
#include "grins/runner.h"

//libMesh
#include "libmesh/exact_solution.h"

int main(int argc, char* argv[])
{
  // Check command line count.
  if( argc < 4 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify libMesh input file, regression file, and regression tolerance." << std::endl;
      exit(1); // TODO: something more sophisticated for parallel runs?
    }

  GRINS::Runner grins(argc,argv);
  grins.init();
  grins.run();

  GRINS::Simulation & sim = grins.get_simulation();

  // Get equation systems to create ExactSolution object
  std::shared_ptr<libMesh::EquationSystems> es = sim.get_equation_system();

  //es->write("foobar.xdr");

  // Create Exact solution object and attach exact solution quantities
  libMesh::ExactSolution exact_sol(*es);

  libMesh::EquationSystems es_ref( es->get_mesh() );

  // Filename of file where comparison solution is stashed
  std::string solution_file = std::string(argv[2]);
  es_ref.read( solution_file );

  exact_sol.attach_reference_solution( &es_ref );

  const GetPot & inputfile = grins.get_input_file();
  std::string system_name = inputfile( "screen-options/system_name", "GRINS" );

  // Compute error and get it in various norms
  exact_sol.compute_error(system_name, "u");
  exact_sol.compute_error(system_name, "v");

  exact_sol.compute_error(system_name, "p");

  double u_l2error = exact_sol.l2_error(system_name, "u");
  double u_h1error = exact_sol.h1_error(system_name, "u");

  double v_l2error = exact_sol.l2_error(system_name, "v");
  double v_h1error = exact_sol.h1_error(system_name, "v");

  double p_l2error = exact_sol.l2_error(system_name, "p");
  double p_h1error = exact_sol.h1_error(system_name, "p");

  int return_flag = 0;

  // This is the tolerance of the iterative linear solver so
  // it's unreasonable to expect anything better than this.
  double tol = atof(argv[3]);

  if( u_l2error > tol || u_h1error > tol ||
      v_l2error > tol || v_h1error > tol ||
      p_l2error > tol || p_h1error > tol  )
    {
      return_flag = 1;

      std::cout << "Tolerance exceeded for thermally driven flow test." << std::endl
                << "tolerance = " << tol << std::endl
                << "u l2 error = " << u_l2error << std::endl
                << "u h1 error = " << u_h1error << std::endl
                << "v l2 error = " << v_l2error << std::endl
                << "v h1 error = " << v_h1error << std::endl
                << "p l2 error = " << p_l2error << std::endl
                << "p h1 error = " << p_h1error << std::endl;
    }

  return return_flag;
}
