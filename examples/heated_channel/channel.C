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
#include "grins/mesh_builder.h"
#include "grins/simulation.h"
#include "parabolic_profile.h"

// GRVY
#ifdef GRINS_HAVE_GRVY
#include "grvy.h"
#endif

// libMesh
#include "libmesh/parallel.h"

class ChannelBCFactory : public GRINS::BoundaryConditionsFactory
{
public:

  ChannelBCFactory( const GetPot& input )
    : GRINS::BoundaryConditionsFactory(input)
  { return; };

  ~ChannelBCFactory(){return;};

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
 
  // MeshBuilder for handling mesh construction
  GRINS::MeshBuilder mesh_builder;

  // PhysicsFactory handles which GRINS::Physics objects to create
  GRINS::PhysicsFactory physics_factory;

  // PhysicsFactory handles which GRINS::Solver to use to solve the problem
  GRINS::SolverFactory solver_factory;

  // VisualizationFactory handles the type of visualization for the simulation
  GRINS::VisualizationFactory vis_factory( libMesh_inputfile );

  ChannelBCFactory bc_factory( libMesh_inputfile );

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
      Real T_init = libMesh_inputfile("Physics/LowMachNavierStokes/T0", 0.0);
      Real p0_init = libMesh_inputfile("Physics/LowMachNavierStokes/p0", 0.0);

      Real& dummy_T  = params.set<Real>("T_init");
      dummy_T = T_init;

      Real& dummy_p0 = params.set<Real>("p0_init");
      dummy_p0 = p0_init;

      system.project_solution( initial_values, NULL, params );
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

Real initial_values( const Point&p, const Parameters &params, 
		     const std::string& , const std::string& unknown_name )
{
  Real value = 0.0;

  if( unknown_name == "T" )
    value = params.get<Real>("T_init");

  else if( unknown_name == "p0" )
    value = params.get<Real>("p0_init");
  
  else if( unknown_name == "u" )
    value = 0.6*p(1)*(1.0-p(1));

  else
    value = 0.0;

  return value;
}

std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > ChannelBCFactory::build_dirichlet( )
{
  GRINS::DBCContainer cont;
  cont.add_var_name( "u" );
  cont.add_bc_id( 1 );
  
  std::tr1::shared_ptr<libMesh::FunctionBase<Number> > vel_func( new GRINS::ParabolicProfile( -0.6, 0.0, 0.0, 0.6, 0.0, 0.0 ) );
    
  cont.set_func( vel_func );


  GRINS::DBCContainer cont2;
  cont2.add_var_name( "v" );
  cont2.add_bc_id( 1 );

  std::tr1::shared_ptr<libMesh::FunctionBase<Number> > vel_func2( new ZeroFunction<Number> );

  cont2.set_func( vel_func2 );

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > mymap;
  
  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::low_mach_navier_stokes,  cont) );

  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::low_mach_navier_stokes,  cont2) );

  return mymap;
}
