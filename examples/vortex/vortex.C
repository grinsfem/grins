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
#include "gaussian_xy_profile.h"

// GRVY
#ifdef HAVE_GRVY
#include "grvy.h"
#endif

// libMesh
#include "parallel.h"

class VortexBCFactory : public GRINS::BoundaryConditionsFactory
{
public:

  VortexBCFactory( const GetPot& input )
    : GRINS::BoundaryConditionsFactory(input)
  { return; };

  ~VortexBCFactory(){return;};

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > build_dirichlet( );
};

// Function for getting initial temperature field
Real initial_values( const Point& p, const Parameters &params, 
		     const std::string& system_name, const std::string& unknown_name );

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
  GRINS::PhysicsFactory physics_factory;

  // PhysicsFactory handles which GRINS::Solver to use to solve the problem
  GRINS::SolverFactory solver_factory( libMesh_inputfile );

  // VisualizationFactory handles the type of visualization for the simulation
  GRINS::VisualizationFactory vis_factory( libMesh_inputfile );

  VortexBCFactory bc_factory( libMesh_inputfile );

  GRINS::Simulation grins( libMesh_inputfile,
			   &physics_factory,
			   &mesh_builder,
			   &solver_factory,
			   &vis_factory,
			   &bc_factory );

  //FIXME: We need to move this to within the Simulation object somehow...
  std::string restart_file = libMesh_inputfile( "restart-options/restart_file", "none" );

  if( restart_file == "none" )
    {
      // Asssign initial temperature value
      std::string system_name = libMesh_inputfile( "screen-options/system_name", "GRINS" );
      std::tr1::shared_ptr<libMesh::EquationSystems> es = grins.get_equation_system();
      const libMesh::System& system = es->get_system(system_name);
      
      Parameters &params = es->parameters;
      Real T_init = libMesh_inputfile("InitialConditions/T0", 0.0);

      Real u_theta = libMesh_inputfile("InitialConditions/u_theta", 0.0);

      Real& dummy_T  = params.set<Real>("T_init");
      dummy_T = T_init;

      Real& dummy_u = params.set<Real>("u_theta");
      dummy_u = u_theta;
      
      std::cout << "==========================================================" << std::endl
		<< " Beginning Solution Projection" << std::endl
		<< "==========================================================" << std::endl;

      system.project_solution( initial_values, NULL, params );

      std::cout << "==========================================================" << std::endl
		<< " Finished Solution Projection" << std::endl
		<< "==========================================================" << std::endl;
    }

#ifdef USE_GRVY_TIMERS
  grvy_timer.EndTimer("Initialize Solver");

  // Attach GRVY timer to solver
  grins.attach_grvy_timer( &grvy_timer );
#endif

  grins.run();

#ifdef USE_GRVY_TIMERS
  grvy_timer.Finalize();
 
  if( Parallel::Communicator_World.rank() == 0 ) grvy_timer.Summarize();
#endif

  return 0;
}

Real initial_values( const Point& p, const Parameters &params, 
		     const std::string& , const std::string& unknown_name )
{
  if( unknown_name == "T" )
    {
      return params.get<Real>("T_init");
    }

  Real r = std::sqrt( p(0)*p(0) + p(1)*p(1) );
  
  Real theta = 0.0;
      
  if( p(0) >= 0.0 )
    theta = std::asin(p(1)/r);
  else
    theta = -std::asin(p(1)/r) + pi;
  
  Real u_theta = params.get<Real>("u_theta");

  Real rc = 0.2;
  Real sigma = 0.2/4;

  Real vert_scaling = 0.0;
  
  if( p(2) >= 0.1 || p(2) <= 0.9 )
    vert_scaling = 1.0;
  
  if( p(2) < 0.1 )
    {
      vert_scaling = 10.0*p(2);
    }

  if( p(2) > 0.9 )
    vert_scaling = 10.0*(1.0 - p(2));

  if( unknown_name == "u" )
    {
      return -vert_scaling*std::sin(theta)*u_theta*std::exp( -(r-rc)*(r-rc)/(2*sigma*sigma));
    }
  else if( unknown_name == "v" )
    {
      return  vert_scaling*std::cos(theta)*u_theta*std::exp( -(r-rc)*(r-rc)/(2*sigma*sigma));
    }

  if( r >= 0.025 )
    {
      if( unknown_name == "w" )
	{
	  return -vert_scaling*0.002506266;
	}
    }
  else
    {
      if( unknown_name == "w" )
	{
	  return vert_scaling*1.0;
	}
    }

  return 0.0;
}

std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > VortexBCFactory::build_dirichlet( )
{
  // Inflow (w-component) at bc-id 1
  GRINS::DBCContainer cont;
  cont.add_var_name( "w" );
  cont.add_bc_id( 1 );
  
  const Real a = 1.0;
  const Real mu = 0.0;
  const Real sigma = 0.025/3.0;
  const Real shift = a*std::exp( -(3.0*sigma - mu)*(3.0*sigma - mu)/(2.0*sigma*sigma) );

  //std::tr1::shared_ptr<libMesh::FunctionBase<Number> > vel_func( new GRINS::GaussianXYProfile( a, mu, sigma, shift ) );
  std::tr1::shared_ptr<libMesh::FunctionBase<Number> > vel_func( new ZeroFunction<Number> );

  cont.set_func( vel_func );

  // Inflow (u,v-components) at bc-id 1
  GRINS::DBCContainer cont2;
  cont2.add_var_name( "u" );
  cont2.add_var_name( "v" );
  cont2.add_bc_id( 1 );

  std::tr1::shared_ptr<libMesh::FunctionBase<Number> > vel_func2( new ZeroFunction<Number> );

  cont2.set_func( vel_func2 );


  // No-flow at bc-id 3. Just setting w-component to 0.
  GRINS::DBCContainer cont3;
  cont3.add_var_name( "w" );
  cont3.add_bc_id( 7 );
  
  // Reuse zero function
  cont3.set_func( vel_func2 );


  // Now pack it all up.
  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > mymap;
  
  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::incompressible_navier_stokes,  cont) );

  //mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::incompressible_navier_stokes,  cont2) );

  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::incompressible_navier_stokes,  cont3) );

  return mymap;
}
