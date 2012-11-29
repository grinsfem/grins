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
#include "simulation.h"
#include "multiphysics_sys.h"

//libMesh
#include "exact_solution.h"

// GRVY
#ifdef HAVE_GRVY
#include "grvy.h"
#endif

Number exact_solution( const Point& p,
		       const Parameters&,   // parameters, not needed
		       const std::string&,  // sys_name, not needed
		       const std::string&); // unk_name, not needed);

Gradient exact_derivative( const Point& p,
			   const Parameters&,   // parameters, not needed
			   const std::string&,  // sys_name, not needed
			   const std::string&); // unk_name, not needed);

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
 
  GRINS::SimulationBuilder sim_builder;

  GRINS::Simulation grins( libMesh_inputfile,
			   sim_builder );

#ifdef USE_GRVY_TIMERS
  grvy_timer.EndTimer("Initialize Solver");

  // Attach GRVY timer to solver
  grins.attach_grvy_timer( &grvy_timer );
#endif

  grins.run();
  
#ifdef USE_GRVY_TIMERS
  grvy_timer.Finalize();
#endif

  std::tr1::shared_ptr<libMesh::EquationSystems> es = grins.get_equation_system();

  // Create Exact solution object and attach exact solution quantities
  ExactSolution exact_sol(*es);

  exact_sol.attach_exact_value(&exact_solution);
  exact_sol.attach_exact_deriv(&exact_derivative);
  
  // Compute error and get it in various norms
  exact_sol.compute_error("GRINS", "u");

  double l2error = exact_sol.l2_error("GRINS", "u");
  double h1error = exact_sol.h1_error("GRINS", "u");

  int return_flag = 0;

  if( l2error > 1.0e-12 || h1error > 1.0e-12 )
    {
      return_flag = 1;

      std::cout << "Tolerance exceeded for Couette flow test." << std::endl
		<< "l2 error = " << l2error << std::endl
		<< "h1 error = " << h1error << std::endl;
    }

 return return_flag;
}

Number exact_solution( const Point& p,
		       const Parameters& params,   // parameters, not needed
		       const std::string&,  // sys_name, not needed
		       const std::string&)  // unk_name, not needed);
{
  const double x = p(0);
  const double y = p(1);
  const double z = p(2);
  
  // Hardcoded to velocity in input file.
  Number f = 10.0*y;

  return f;
}

Gradient exact_derivative( const Point& p,
			   const Parameters& params,   // parameters, not needed
			   const std::string&,  // sys_name, not needed
			   const std::string&)  // unk_name, not needed);
{
  const double x = p(0);
  const double y = p(1);
  const double z = p(2);

  Gradient g;

  // Hardcoded to velocity in input file.
  g(0) = 0.0;
  g(1) = 10.0;
  g(2) = 0.0;

  return g;
}
