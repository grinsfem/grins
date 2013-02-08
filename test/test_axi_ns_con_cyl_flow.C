//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "grins_config.h"

#include <iostream>

// GRINS
#include "grins/mesh_builder.h"
#include "grins/simulation.h"
#include "grins/simulation_builder.h" 
#include "grins/multiphysics_sys.h"
#include "grins/concentric_cylinder_profile.h"

//libMesh
#include "libmesh/exact_solution.h"
#include "libmesh/zero_function.h"

// GRVY
#ifdef GRINS_HAVE_GRVY
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

class AxiConCylBCFactory : public GRINS::BoundaryConditionsFactory
{
public:

  AxiConCylBCFactory( )
    : GRINS::BoundaryConditionsFactory()
  { return; };

  ~AxiConCylBCFactory(){return;};

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > build_dirichlet( );
};

int main(int argc, char* argv[]) 
{
#ifdef GRINS_USE_GRVY_TIMERS
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

#ifdef GRINS_USE_GRVY_TIMERS
  grvy_timer.BeginTimer("Initialize Solver");
#endif

  // Initialize libMesh library.
  LibMeshInit libmesh_init(argc, argv);
 
  GRINS::SimulationBuilder sim_builder;

  std::tr1::shared_ptr<AxiConCylBCFactory> bc_factory( new AxiConCylBCFactory );

  sim_builder.attach_bc_factory(bc_factory);

  GRINS::Simulation grins( libMesh_inputfile,
			   sim_builder );

#ifdef GRINS_USE_GRVY_TIMERS
  grvy_timer.EndTimer("Initialize Solver");

  // Attach GRVY timer to solver
  grins.attach_grvy_timer( &grvy_timer );
#endif

  // Do solve here
  grins.run();

  // Get equation systems to create ExactSolution object
  std::tr1::shared_ptr<EquationSystems> es = grins.get_equation_system();

  // Create Exact solution object and attach exact solution quantities
  ExactSolution exact_sol(*es);

  exact_sol.attach_exact_value(&exact_solution);
  exact_sol.attach_exact_deriv(&exact_derivative);
  
  // Compute error and get it in various norms
  exact_sol.compute_error("GRINS", "z_vel");

  double l2error = exact_sol.l2_error("GRINS", "z_vel");
  double h1error = exact_sol.h1_error("GRINS", "z_vel");
  
  int return_flag = 0;

  if( l2error > 1.0e-10 || h1error > 4.0e-7 )
    {
      return_flag = 1;

      std::cout << "Tolerance exceeded for velocity in Poiseuille test." << std::endl
		<< "l2 error = " << l2error << std::endl
		<< "h1 error = " << h1error << std::endl;
    }

  /*
  // Compute error and get it in various norms
  exact_sol.compute_error("GRINS", "p");

  l2error = exact_sol.l2_error("GRINS", "p");
  h1error = exact_sol.h1_error("GRINS", "p");

  if( l2error > 1.0e-11 || h1error > 1.0e-11 )
    {
      return_flag = 1;

      std::cout << "Tolerance exceeded for pressure in Poiseuille test." << std::endl
		<< "l2 error = " << l2error << std::endl
		<< "h1 error = " << h1error << std::endl;
    }
  */

  return return_flag;
}

std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > AxiConCylBCFactory::build_dirichlet( )
{
  GRINS::DBCContainer cont;
  cont.add_var_name( "z_vel" );
  cont.add_bc_id( 0 );
  cont.add_bc_id( 2 );
  
  std::tr1::shared_ptr<libMesh::FunctionBase<Number> > vel_func( new GRINS::ConcentricCylinderProfile );
    
  cont.set_func( vel_func );


  GRINS::DBCContainer cont2;
  cont2.add_var_name( "z_vel" );
  cont2.add_bc_id( 0 );
  cont2.add_bc_id( 2 );

  std::tr1::shared_ptr<libMesh::FunctionBase<Number> > vel_func2( new ZeroFunction<Number> );

  cont2.set_func( vel_func2 );

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > mymap;
  
  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::axisymmetric_incomp_navier_stokes,  cont) );

  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::axisymmetric_incomp_navier_stokes,  cont2) );

  return mymap;
}

Number exact_solution( const Point& p,
		       const Parameters& /*params*/,   // parameters, not needed
		       const std::string& /*sys_name*/,  // sys_name, not needed
		       const std::string& var )  // unk_name, not needed);
{
  const double r = p(0);
  
  const double r0 = 1.0;
  const double r1 = 2.0;
  const double u0 = 2.0;

  Number f;
  // Hardcoded to velocity in input file.
  if( var == "z_vel" ) f = u0*std::log( r1/r )/std::log( r1/r0 );

  return f;
}

Gradient exact_derivative( const Point& p,
			   const Parameters& /*params*/,   // parameters, not needed
			   const std::string& /*sys_name*/,  // sys_name, not needed
			   const std::string& var)  // unk_name, not needed);
{
  const double r = p(0);

  const double r0 = 1.0;
  const double r1 = 2.0;
  const double u0 = 2.0;

  Gradient g;

  // Hardcoded to velocity in input file.
  if( var == "z_vel" )
    {
      g(0) = -u0/std::log(r1/r0)*r/r1*r1/(r*r);
      g(1) = 0.0;
    }

  return g;
}
