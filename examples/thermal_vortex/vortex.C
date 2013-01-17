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
#include "grins/simulation.h"
#include "grins/gaussian_xy_profile.h"

// GRVY
#ifdef GRINS_HAVE_GRVY
#include "grvy.h"
#endif

// libMesh
#include "libmesh/parallel.h"

class VortexBCFactory : public GRINS::BoundaryConditionsFactory
{
public:

  VortexBCFactory( )
    : GRINS::BoundaryConditionsFactory()
  { return; };

  ~VortexBCFactory(){return;};

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > build_dirichlet( );
};

// Function for getting initial temperature field
Real initial_values( const Point& p, const Parameters &params, 
		     const std::string& system_name, const std::string& unknown_name );

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

  std::tr1::shared_ptr<VortexBCFactory> bc_factory( new VortexBCFactory );

  sim_builder.attach_bc_factory( bc_factory );

  GRINS::Simulation grins( libMesh_inputfile,
			   sim_builder );

  //FIXME: We need to move this to within the Simulation object somehow...
  std::string restart_file = libMesh_inputfile( "restart-options/restart_file", "none" );

  if( restart_file == "none" )
    {
      // Asssign initial temperature value
      std::string system_name = libMesh_inputfile( "screen-options/system_name", "GRINS" );
      std::tr1::shared_ptr<libMesh::EquationSystems> es = grins.get_equation_system();
      const libMesh::System& system = es->get_system(system_name);
      
      Parameters &params = es->parameters;
      Real T_base = libMesh_inputfile("InitialConditions/T_base", 0.0);
      Real T_top = libMesh_inputfile("InitialConditions/T_top", 0.0);

      Real u_theta = libMesh_inputfile("InitialConditions/u_theta", 0.0);

      Real jet_width = libMesh_inputfile("InitialConditions/jet_width", 0.0);

      Real w0 = libMesh_inputfile("InitialConditions/w0", 0.0);

      Real& dummy_Tb  = params.set<Real>("T_base");
      dummy_Tb = T_base;

      Real& dummy_Tt  = params.set<Real>("T_top");
      dummy_Tt = T_top;

      Real& dummy_u = params.set<Real>("u_theta");
      dummy_u = u_theta;

      Real& dummy_jw = params.set<Real>("jet_width");
      dummy_jw = jet_width;
      
      Real& dummy_w0 = params.set<Real>("w0");
      dummy_w0 = w0;

      std::cout << "==========================================================" << std::endl
		<< " Beginning Solution Projection" << std::endl
		<< "==========================================================" << std::endl;

      system.project_solution( initial_values, NULL, params );

      std::cout << "==========================================================" << std::endl
		<< " Finished Solution Projection" << std::endl
		<< "==========================================================" << std::endl;
    }

#ifdef GRINS_USE_GRVY_TIMERS
  grvy_timer.EndTimer("Initialize Solver");

  // Attach GRVY timer to solver
  grins.attach_grvy_timer( &grvy_timer );
#endif

  grins.run();

#ifdef GRINS_USE_GRVY_TIMERS
  grvy_timer.Finalize();
 
  if( Parallel::Communicator_World.rank() == 0 ) grvy_timer.Summarize();
#endif

  return 0;
}

Real initial_values( const Point& p, const Parameters &params, 
		     const std::string& , const std::string& unknown_name )
{
  Real x = p(0);
  Real y = p(1);
  Real z = p(2);

  std::string uvar = "u";
  std::string wvar = "w";

  /*
  x = p(2);
  z = p(0);
  
  uvar = "w";
  wvar = "u";
  */
  

  if( unknown_name == "T" )
    {
      const Real T_base = params.get<Real>("T_base");
      const Real T_top = params.get<Real>("T_top");
      return (T_top - T_base)*z + T_base;
    }

  Real r = std::sqrt( x*x + y*y );
  
  Real theta = 0.0;
      
  if( r > 1.0e-6)
    {
      if( x >= 0.0 )
	theta = std::asin(y/r);
      else
	theta = -std::asin(y/r) + pi;
    }

  Real u_theta = params.get<Real>("u_theta");

  Real rc = 0.2;
  Real sigma = 0.2/4;

  Real vert_scaling = 0.0;
  
  if( z >= 0.1 || z <= 0.9 )
    vert_scaling = 1.0;
  
  if( z < 0.1 )
    {
      vert_scaling = 10.0*z;
    }

  if( z > 0.9 )
    vert_scaling = 10.0*(1.0 - z);

  if( unknown_name == uvar )
    {
      return -vert_scaling*std::sin(theta)*u_theta*std::exp( -(r-rc)*(r-rc)/(2*sigma*sigma));
    }
  else if( unknown_name == "v" )
    {
      return  vert_scaling*std::cos(theta)*u_theta*std::exp( -(r-rc)*(r-rc)/(2*sigma*sigma));
    }

  const Real jet_width = params.get<Real>("jet_width");
  const Real w0 = params.get<Real>("w0");

  if( std::fabs(x) <= jet_width/2.0 && std::fabs(y) <= jet_width/2.0 )
    {
      if( unknown_name == wvar )
	{
	  return vert_scaling*w0;
	}
    }
  else
    {
      if( unknown_name == wvar )
	{
	  const Real A0 = jet_width*jet_width;
	  const Real A1 = 1.0-A0;
	  return -vert_scaling*w0*A0/A1;
	}
    }

  return 0.0;
}

std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > VortexBCFactory::build_dirichlet( )
{
  std::tr1::shared_ptr<libMesh::FunctionBase<Number> > vel_func( new ZeroFunction<Number> );

  // No-flow at top.
  GRINS::DBCContainer cont3;
  cont3.add_var_name( "w" );
  cont3.add_bc_id( 5 );
  //cont3.add_bc_id( 0 );
  
  // Reuse zero function
  cont3.set_func( vel_func );


  // No-flow at xmin,xmax face
  GRINS::DBCContainer cont4;
  cont4.add_var_name( "u" );
  cont4.add_bc_id( 2 );
  cont4.add_bc_id( 4 );
  // Reuse zero function
  cont4.set_func( vel_func );

  // No-flow at ymin,ymax face
  GRINS::DBCContainer cont5;
  cont5.add_var_name( "v" );
  cont5.add_bc_id( 3 );
  cont5.add_bc_id( 1 );
  // Reuse zero function
  cont5.set_func( vel_func );

  // Now pack it all up.
  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > mymap;
  
  //mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::incompressible_navier_stokes,  cont) );

  //mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::incompressible_navier_stokes,  cont2) );

  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::incompressible_navier_stokes,  cont3) );
  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::incompressible_navier_stokes,  cont4) );
  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::incompressible_navier_stokes,  cont5) );

  return mymap;
}
