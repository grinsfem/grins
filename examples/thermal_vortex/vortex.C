//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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
#include "grins/bc_factory.h"
#include "grins/gaussian_xy_profile.h"
#include "grins/simulation.h"
#include "grins/simulation_builder.h"

// GRVY
#ifdef GRINS_HAVE_GRVY
#include "grvy.h"
#endif

// libMesh
#include "libmesh/parallel.h"
#include "libmesh/zero_function.h"

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
libMesh::Real
initial_values( const libMesh::Point& p, const libMesh::Parameters &params, 
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
  libMesh::LibMeshInit libmesh_init(argc, argv);
 
  GRINS::SimulationBuilder sim_builder;

  std::tr1::shared_ptr<VortexBCFactory> bc_factory( new VortexBCFactory );

  sim_builder.attach_bc_factory( bc_factory );

  GRINS::Simulation grins( libMesh_inputfile,
			   sim_builder,
                           libmesh_init.comm() );

  //FIXME: We need to move this to within the Simulation object somehow...
  std::string restart_file = libMesh_inputfile( "restart-options/restart_file", "none" );

  if( restart_file == "none" )
    {
      // Asssign initial temperature value
      std::string system_name = libMesh_inputfile( "screen-options/system_name", "GRINS" );
      std::tr1::shared_ptr<libMesh::EquationSystems> es = grins.get_equation_system();
      const libMesh::System& system = es->get_system(system_name);
      
      libMesh::Parameters &params = es->parameters;
      libMesh::Real T_base = libMesh_inputfile("InitialConditions/T_base", 0.0);
      libMesh::Real T_top = libMesh_inputfile("InitialConditions/T_top", 0.0);

      libMesh::Real u_theta = libMesh_inputfile("InitialConditions/u_theta", 0.0);

      libMesh::Real jet_width = libMesh_inputfile("InitialConditions/jet_width", 0.0);

      libMesh::Real w0 = libMesh_inputfile("InitialConditions/w0", 0.0);

      libMesh::Real& dummy_Tb  = params.set<libMesh::Real>("T_base");
      dummy_Tb = T_base;

      libMesh::Real& dummy_Tt  = params.set<libMesh::Real>("T_top");
      dummy_Tt = T_top;

      libMesh::Real& dummy_u = params.set<libMesh::Real>("u_theta");
      dummy_u = u_theta;

      libMesh::Real& dummy_jw = params.set<libMesh::Real>("jet_width");
      dummy_jw = jet_width;
      
      libMesh::Real& dummy_w0 = params.set<libMesh::Real>("w0");
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

libMesh::Real
initial_values( const libMesh::Point& p, const libMesh::Parameters &params, 
		const std::string& , const std::string& unknown_name )
{
  libMesh::Real x = p(0);
  libMesh::Real y = p(1);
  libMesh::Real z = p(2);

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
      const libMesh::Real T_base = params.get<libMesh::Real>("T_base");
      const libMesh::Real T_top = params.get<libMesh::Real>("T_top");
      return (T_top - T_base)*z + T_base;
    }

  libMesh::Real r = std::sqrt( x*x + y*y );
  
  libMesh::Real theta = 0.0;
      
  if( r > 1.0e-6)
    {
      if( x >= 0.0 )
	theta = std::asin(y/r);
      else
	theta = -std::asin(y/r) + libMesh::pi;
    }

  libMesh::Real u_theta = params.get<libMesh::Real>("u_theta");

  libMesh::Real rc = 0.2;
  libMesh::Real sigma = 0.2/4;

  libMesh::Real vert_scaling = 0.0;
  
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

  const libMesh::Real jet_width = params.get<libMesh::Real>("jet_width");
  const libMesh::Real w0 = params.get<libMesh::Real>("w0");

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
	  const libMesh::Real A0 = jet_width*jet_width;
	  const libMesh::Real A1 = 1.0-A0;
	  return -vert_scaling*w0*A0/A1;
	}
    }

  return 0.0;
}

std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > VortexBCFactory::build_dirichlet( )
{
  std::tr1::shared_ptr<libMesh::FunctionBase<libMesh::Number> >
    vel_func( new libMesh::ZeroFunction<libMesh::Number> );

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
