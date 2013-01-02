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
// $Id: injection.C 35065 2012-12-05 04:46:12Z pbauman $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#include "grins_config.h"

#include <iostream>

// GRINS
#include "simulation.h"
#include "constant_with_exp_layer.h"

// GRVY
#ifdef GRINS_HAVE_GRVY
#include "grvy.h"
#endif

// libMesh
#include "parallel.h"

class BunsenBCFactory : public GRINS::BoundaryConditionsFactory
{
public:

  BunsenBCFactory( )
    : GRINS::BoundaryConditionsFactory()
  { return; };

  ~BunsenBCFactory(){return;};

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

  std::tr1::shared_ptr<BunsenBCFactory> bc_factory( new BunsenBCFactory );

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
      Real T_init = libMesh_inputfile("InitialConditions/T0", 0.0);

      Real& dummy_T  = params.set<Real>("T_init");
      dummy_T = T_init;

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

Real initial_values( const Point& /*p*/, const Parameters &params, 
		     const std::string& , const std::string& unknown_name )
{
  Real value = 0.0;

  if( unknown_name == "T" )
    {
      value = params.get<Real>("T_init");
    }
  else
    {
      value = 0.0;
    }

  return value;
}

std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > BunsenBCFactory::build_dirichlet( )
{
  const Real delta = 0.0005;
  const Real u0 = 1.2;

  GRINS::DBCContainer cont;
  {
    cont.add_var_name( "v" );
    cont.add_bc_id( 1 );
    
    const Real factor = 1.0;
    const Real r0 = 0.002;

    std::tr1::shared_ptr<libMesh::FunctionBase<Number> > vel_func( new GRINS::ConstantWithExponentialLayer( u0, factor, r0, delta ) );
    
    cont.set_func( vel_func );
  }

  GRINS::DBCContainer cont2;
  {
    cont2.add_var_name( "u" );
    cont2.add_bc_id( 1 );
    cont2.add_bc_id( 3 );
    
    std::tr1::shared_ptr<libMesh::FunctionBase<Number> > vel_func( new ZeroFunction<Number> );

    cont2.set_func( vel_func );
  }

  GRINS::DBCContainer cont3;
  {
    cont3.add_var_name( "v" );
    cont3.add_bc_id( 3 );
    
    const Real factor = -1.0;
    const Real r0 = 0.0025;

    std::tr1::shared_ptr<libMesh::FunctionBase<Number> > vel_func( new GRINS::ConstantWithExponentialLayer( u0, factor, r0, delta ) );
    
    cont3.set_func( vel_func );
  }

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > mymap;
  
  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::reacting_low_mach_navier_stokes,  cont) );

  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::reacting_low_mach_navier_stokes,  cont2) );

  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::reacting_low_mach_navier_stokes,  cont3) );

  return mymap;
}
