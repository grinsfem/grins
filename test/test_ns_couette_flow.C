//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "config.h"

#include <iostream>

// GRINS stuff
#include "grins_mesh_manager.h"
#include "grins_solver.h"

// System types that we might want to instantiate
#include "multiphysics_sys.h"

#include "exact_solution.h"

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

  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify libMesh input file." << std::endl;
      exit(1); // TODO: something more sophisticated for parallel runs?
    }

  // Initialize libMesh library.
  LibMeshInit libmesh_init(argc, argv);
  
  // Create mesh manager object.
  GRINS::MeshManager meshmanager;
  
  // Create solver object.
  GRINS::Solver<GRINS::MultiphysicsSystem> solver;

  { // Artificial block to destroy objects associated with reading the input once we've read it in.

    // libMesh input file should be first argument
    std::string libMesh_input_filename = argv[1];
    
    // Create our GetPot object. TODO: Finalize decision of GRVY vs. GetPot input.
    GetPot libMesh_inputfile( libMesh_input_filename );
    
    // Read solver options
    solver.read_input_options( libMesh_inputfile );

    // Read mesh options
    meshmanager.read_input_options( libMesh_inputfile );

    // Setup and initialize system so system can read it's relavent options
    meshmanager.build_mesh();
    solver.set_mesh( meshmanager.get_mesh() );
    solver.initialize_system( "GRINS", libMesh_inputfile );

  } //Should be done reading input, so we kill the GetPot object.

  // Do solve here
  solver.solve();

  // Get equation systems to create ExactSolution object
  GRINS::MultiphysicsSystem* system = solver.get_system();

  EquationSystems & es = system->get_equation_systems ();

  // Create Exact solution object and attach exact solution quantities
  ExactSolution exact_sol(es);

  exact_sol.attach_exact_value(&exact_solution);
  exact_sol.attach_exact_deriv(&exact_derivative);
  
  // Compute error and get it in various norms
  exact_sol.compute_error("GRINS", "u");

  double l2error = exact_sol.l2_error("GRINS", "u");
  double h1error = exact_sol.h1_error("GRINS", "u");

  int return_flag = 0;

  if( l2error > 6.0e-10 || h1error > 3.0e-9 )
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
