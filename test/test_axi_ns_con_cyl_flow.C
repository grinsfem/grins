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

#include "config.h"

#include <iostream>

// GRINS
#include "mesh_builder.h"
#include "simulation.h"
#include "multiphysics_sys.h"
#include "concentric_cylinder_profile.h"

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

class MyBCFactory : public GRINS::BoundaryConditionsFactory
{
public:

  MyBCFactory( const GetPot& input )
    : GRINS::BoundaryConditionsFactory(input)
  { return; };

  ~MyBCFactory(){return;};

  std::map< std::string, GRINS::DBCContainer >
  build_dirichlet( libMesh::EquationSystems& equation_system );
};

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
 
  // MeshBuilder for handling mesh construction
  GRINS::MeshBuilder mesh_builder( libMesh_inputfile );

  // PhysicsFactory handles which GRINS::Physics objects to create
  GRINS::PhysicsFactory physics_factory( libMesh_inputfile );

  // PhysicsFactory handles which GRINS::Solver to use to solve the problem
  GRINS::SolverFactory solver_factory( libMesh_inputfile );

  // VisualizationFactory handles the type of visualization for the simulation
  GRINS::VisualizationFactory vis_factory( libMesh_inputfile );

  MyBCFactory bc_factory( libMesh_inputfile );

  GRINS::Simulation grins( libMesh_inputfile,
			   &physics_factory,
			   &mesh_builder,
			   &solver_factory,
			   &vis_factory,
			   &bc_factory );

#ifdef USE_GRVY_TIMERS
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

std::map< std::string, GRINS::DBCContainer > MyBCFactory::build_dirichlet( libMesh::EquationSystems& es )
{
  const libMesh::System& system = es.get_system("GRINS");
  GRINS::VariableIndex z_var = system.variable_number("z_vel");

  std::tr1::shared_ptr<GRINS::DirichletFuncObj> inflow( new GRINS::ConcentricCylinderProfile(z_var) );

  GRINS::DirichletBCsMap dbc_map;
  dbc_map.insert( GRINS::DBCMapPair( z_var,inflow ) );

  GRINS::DBCContainer dbc_container;
  dbc_container.insert( GRINS::DBCContainerPair( 0, dbc_map ) );
  dbc_container.insert( GRINS::DBCContainerPair( 2, dbc_map ) );

  std::map< std::string, GRINS::DBCContainer > dbcs;

  dbcs.insert( std::pair< std::string, GRINS::DBCContainer >( "AxisymmetricIncompressibleNavierStokes", dbc_container ) );

  return dbcs;
}

Number exact_solution( const Point& p,
		       const Parameters& params,   // parameters, not needed
		       const std::string& sys_name,  // sys_name, not needed
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
			   const Parameters& params,   // parameters, not needed
			   const std::string& sys_name,  // sys_name, not needed
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
