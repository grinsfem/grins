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
// $Id: injection.C 35065 2012-12-05 04:46:12Z pbauman $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#include "grins_config.h"

#include <iostream>

// Bunsen
#include "constant_with_exp_layer.h"
#include "bunsen_source.h"
#include "ignite_initial_guess.h"

// GRINS
#include "grins/simulation.h"
#include "grins/simulation_builder.h"
#include "grins/physics_factory.h"
#include "grins/cantera_singleton.h"
#include "grins/bc_factory.h"
#include "grins/physics_factory.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/zero_function.h"
#include "libmesh/parallel.h"

// Cantera
#ifdef GRINS_HAVE_CANTERA
#include "cantera/equilibrium.h"
#endif // GRINS_HAVE_CANTERA

// GRVY
#ifdef GRINS_HAVE_GRVY
#include "grvy.h"
#endif

class BunsenBCFactory : public GRINS::BoundaryConditionsFactory
{
public:

  BunsenBCFactory( )
    : GRINS::BoundaryConditionsFactory()
  { return; };

  ~BunsenBCFactory(){return;};

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > build_dirichlet( );
};

class BunsenPhysicsFactory : public GRINS::PhysicsFactory
{
public:

  BunsenPhysicsFactory() : GRINS::PhysicsFactory(){};
  virtual ~BunsenPhysicsFactory(){};

protected:
  
  virtual void add_physics( const GetPot& input,
			    const std::string& physics_to_add,
			    GRINS::PhysicsList& physics_list );

  virtual void check_physics_consistency( const GRINS::PhysicsList& physics_list );

};

// Function for getting initial temperature field
libMesh::Real initial_values( const Point& p, const Parameters &params, 
		     const std::string& system_name, const std::string& unknown_name );

static libMesh::Point p_old = libMesh::Point( -1000000000.0, -1000000000.0, -1000000000.0 );

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

  std::tr1::shared_ptr<GRINS::PhysicsFactory> physics_factory( new BunsenPhysicsFactory );
  sim_builder.attach_physics_factory( physics_factory );

  GRINS::Simulation grins( libMesh_inputfile,
			   sim_builder );

  //FIXME: We need to move this to within the Simulation object somehow...
  std::string restart_file = libMesh_inputfile( "restart-options/restart_file", "none" );

  // If we are "cold starting", setup the flow field.
  if( restart_file == "none" )
    {
      // Asssign initial temperature value
      std::string system_name = libMesh_inputfile( "screen-options/system_name", "GRINS" );
      std::tr1::shared_ptr<libMesh::EquationSystems> es = grins.get_equation_system();
      const libMesh::System& system = es->get_system(system_name);
      
      Parameters &params = es->parameters;

      libMesh::Real& T_init = params.set<libMesh::Real>("T_init");
      T_init = libMesh_inputfile("InitialConditions/T0", 0.0);

      libMesh::Real& p0 = params.set<libMesh::Real>("p0");
      p0 = libMesh_inputfile("Physics/ReactingLowMachNavierStokes/p0", 1.0e5);

#ifdef GRINS_HAVE_CANTERA
      Cantera::IdealGasMix& cantera = GRINS::CanteraSingleton::cantera_instance(libMesh_inputfile);
#endif

      libMesh::Real& w_H2 = params.set<libMesh::Real>( "w_H2" );
      w_H2 = libMesh_inputfile( "Physics/ReactingLowMachNavierStokes/bound_species_1", 0.0, 0 );

      libMesh::Real& w_O2 = params.set<libMesh::Real>( "w_O2" );
      w_O2 = libMesh_inputfile( "Physics/ReactingLowMachNavierStokes/bound_species_1", 0.0, 1 );

      libMesh::Real& w_H2O = params.set<libMesh::Real>( "w_H2O" );
      w_H2O = libMesh_inputfile( "Physics/ReactingLowMachNavierStokes/bound_species_1", 0.0, 2 );

      libMesh::Real& w_H = params.set<libMesh::Real>( "w_H" );
      w_H = libMesh_inputfile( "Physics/ReactingLowMachNavierStokes/bound_species_1", 0.0, 3 );

      libMesh::Real& w_O = params.set<libMesh::Real>( "w_O" );
      w_O = libMesh_inputfile( "Physics/ReactingLowMachNavierStokes/bound_species_1", 0.0, 4 );

      libMesh::Real& w_OH = params.set<libMesh::Real>( "w_OH" );
      w_OH = libMesh_inputfile( "Physics/ReactingLowMachNavierStokes/bound_species_1", 0.0, 5 );

      libMesh::Real& w_HO2 = params.set<libMesh::Real>( "w_HO2" );
      w_HO2 = libMesh_inputfile( "Physics/ReactingLowMachNavierStokes/bound_species_1", 0.0, 6 );

      libMesh::Real& w_H2O2 = params.set<libMesh::Real>( "w_H2O2" );
      w_H2O2 = libMesh_inputfile( "Physics/ReactingLowMachNavierStokes/bound_species_1", 0.0, 7 );

      libMesh::Real& w_N2 = params.set<libMesh::Real>( "w_N2" );
      w_N2 = libMesh_inputfile( "Physics/ReactingLowMachNavierStokes/bound_species_1", 0.0, 8 );

      std::cout << "==============================================" << std::endl;
      std::cout << "Projecting Solution." << std::endl;
      std::cout << "==============================================" << std::endl;
      system.project_solution( initial_values, NULL, params );
      std::cout << "==============================================" << std::endl;
      std::cout << "Done Projecting Solution!" << std::endl;
      std::cout << "==============================================" << std::endl;
    }

  /* If we're restarting to try and get ignition, then we need to setup a
     "restart" system and using the IgniteInitalGuess functor to do the projection
     on the "real" system. */
  if( libMesh_inputfile( "restart-options/ignition", false ) && 
      restart_file != std::string("none") )
    {
      std::string system_name = libMesh_inputfile( "screen-options/system_name", "GRINS" );

      /*
      GetPot restart_input( "bunsen_restart.in" );
      
      GRINS::Simulation restart_sim( restart_input,
				     sim_builder );
      std::tr1::shared_ptr<libMesh::EquationSystems> restart_es = restart_sim.get_equation_system();
      libMesh::System& restart_system = restart_es->get_system(system_name);
      GRINS::MultiphysicsSystem& restart_ms_system = libmesh_cast_ref<GRINS::MultiphysicsSystem&>( restart_system );
      */

      std::tr1::shared_ptr<libMesh::EquationSystems> es = grins.get_equation_system();
      libMesh::System& system = es->get_system(system_name);
      GRINS::MultiphysicsSystem& ms_system = libmesh_cast_ref<GRINS::MultiphysicsSystem&>( system );

      Bunsen::IgniteInitialGuess<libMesh::Real> ignite( libMesh_inputfile, ms_system, 
					       ms_system );

      es->reinit();
      
      std::cout << "==============================================" << std::endl;
      std::cout << "Projecting Solution." << std::endl;
      std::cout << "==============================================" << std::endl;
      ms_system.project_solution( &ignite );
      std::cout << "==============================================" << std::endl;
      std::cout << "Done Projecting Solution!" << std::endl;
      std::cout << "==============================================" << std::endl;
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

libMesh::Real initial_values( const Point& p, const Parameters &params, 
		     const std::string& , const std::string& unknown_name )
{
  libMesh::Real value = 0.0;

  const libMesh::Real r = p(0);
  const libMesh::Real z = p(1);
  libMesh::Real T = 0.0;

  /*
  if( z > 0.02 && z <= 0.04 )
    T = (298 - params.get<libMesh::Real>("T_init"))/0.02*(z-0.02) + params.get<libMesh::Real>("T_init");
  else if( 0.005 <= z && z <= 0.02 )
    T = params.get<libMesh::Real>("T_init");
  else
  */
    T = 298.0;

  libMesh::Real p0 = 0.0;
  p0 = params.get<libMesh::Real>("p0");

  /*
  if( unknown_name.find( "w_" ) != std::string::npos )
    {
      Cantera::IdealGasMix& cantera = GRINS::CanteraSingleton::cantera_instance();
      cantera.setState_TP(T,p0);
      Cantera::equilibrate( cantera, "TP" );
    }
  */

  if( unknown_name == "T" )
    {
      value = T;
    }
  else if( unknown_name == "w_H2" ||
	   unknown_name == "w_O2" ||
	   unknown_name == "w_H2O" || 
	   unknown_name == "w_H" ||
	   unknown_name == "w_O" ||
	   unknown_name == "w_OH" ||
	   unknown_name == "w_HO2" ||
	   unknown_name == "w_H2O2"||
	   unknown_name == "w_N2" )
    {
      //if( T == 1000.0 )
	{
	  if( unknown_name == "w_N2" )
	    {
	      value = 0.8;
	    }
	  else if( unknown_name == "w_O2" )
	    {
	      value = 0.2;
	    }
	  else if( unknown_name.find("w_") != std::string::npos )
	    value = 0.0;     
	}
      /*
      else
	{
	  std::vector<libMesh::Real> Y(9,0.0);
	  Cantera::IdealGasMix& cantera = GRINS::CanteraSingleton::cantera_instance();
	  cantera.getMassFractions(&Y[0]);
	  value = Y[ cantera.speciesIndex( unknown_name ) ];
	  
	  std::cout << "T = " << T << std::endl;
	  for( unsigned int s = 0; s < 9; s++ )
	    {
	      std::cout << "Y[" << s << "] = " << Y[s] << std::endl;
	    }
	  
	}
      */
    }
  else
    {
      value = 0.0;
    }

  return value;
}

std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > BunsenBCFactory::build_dirichlet( )
{
  const libMesh::Real delta = 0.0005;
  const libMesh::Real u0 = 0.12;

  GRINS::DBCContainer cont;
  {
    cont.add_var_name( "v" );
    cont.add_bc_id( 1 );
    
    const libMesh::Real factor = 1.0;
    const libMesh::Real r0 = 0.002;

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
    
    const libMesh::Real factor = -1.0;
    const libMesh::Real r0 = 0.0025;

    std::tr1::shared_ptr<libMesh::FunctionBase<Number> > vel_func( new GRINS::ConstantWithExponentialLayer( u0, factor, r0, delta ) );
    
    cont3.set_func( vel_func );
  }

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > mymap;
  
  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::reacting_low_mach_navier_stokes,  cont) );

  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::reacting_low_mach_navier_stokes,  cont2) );

  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::reacting_low_mach_navier_stokes,  cont3) );

  return mymap;
}

void BunsenPhysicsFactory::add_physics( const GetPot& input,
					const std::string& physics_to_add,
					GRINS::PhysicsList& physics_list )
{
  if( physics_to_add == "BunsenSource" )
    {
      physics_list[physics_to_add] = 
	std::tr1::shared_ptr<GRINS::Physics>( new Bunsen::BunsenSource(physics_to_add,input) );
    }
  else
    {
      GRINS::PhysicsFactory::add_physics(input, physics_to_add, physics_list );
    }

  return;
}

void BunsenPhysicsFactory::check_physics_consistency( const GRINS::PhysicsList& physics_list )
{
  for( GRINS::PhysicsListIter physics = physics_list.begin();
       physics != physics_list.end();
       physics++ )
    {
      if( physics->first == "BunsenSource" )
	{
	  if( physics_list.find(GRINS::reacting_low_mach_navier_stokes) == physics_list.end() )
	    {
	      this->physics_consistency_error( physics->first, GRINS::reacting_low_mach_navier_stokes  );
	    }
	}
    }
  
  GRINS::PhysicsFactory::check_physics_consistency( physics_list );

  return;
}
