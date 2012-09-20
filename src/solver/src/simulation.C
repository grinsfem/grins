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

#include "simulation.h"

GRINS::Simulation::Simulation( const GetPot& input,
			       GRINS::PhysicsFactory* physics_factory,
			       GRINS::MeshBuilder* mesh_builder,
			       GRINS::SolverFactory* solver_factory,
			       GRINS::VisualizationFactory* vis_factory,
			       GRINS::BoundaryConditionsFactory* bc_factory )
  :  _mesh( mesh_builder->build() ),
     _equation_system( new libMesh::EquationSystems( *_mesh ) ),
     _solver( solver_factory->build() ),
     _system_name( input("screen-options/system_name", "GRINS" ) ),
     _multiphysics_system( &(_equation_system->add_system<GRINS::MultiphysicsSystem>( _system_name )) ),
     _vis( vis_factory->build() ),
     _print_mesh_info( input("screen-options/print_mesh_info", false ) ),
     _print_log_info( input("screen-options/print_log_info", false ) ),
     _print_equation_system_info( input("screen-options/print_equation_system_info", false ) ),
     _output_vis( input("vis-options/output_vis", false ) ),
     _output_residual( input( "vis-options/output_residual", false ) )
{
  // Only print libMesh logging info if the user requests it
  libMesh::perflog.disable_logging();
  if( this->_print_log_info ) libMesh::perflog.enable_logging();

  GRINS::PhysicsList physics_list = physics_factory->build(input);

  _multiphysics_system->attach_physics_list( physics_list );

  _multiphysics_system->read_input_options( input );

  // This *must* be done before equation_system->init
  if( bc_factory )
    {
      this->attach_dirichlet_bc_funcs( bc_factory->build_dirichlet(), _multiphysics_system );
    }

  _solver->initialize( input, _equation_system, _multiphysics_system );

  // This *must* be done after equation_system->init
  if( bc_factory )
    {
      this->attach_neumann_bc_funcs( bc_factory->build_neumann( *_equation_system ), _multiphysics_system );
    }

  this->check_for_restart( input );

  return;
}

GRINS::Simulation::~Simulation()
{
  return;
}

void GRINS::Simulation::run()
{
  this->print_sim_info();
  
  _solver->solve(  _multiphysics_system, _equation_system, _vis, _output_vis, _output_residual );

  return;
}

void GRINS::Simulation::print_sim_info()
{
  // Print mesh info if the user wants it
  if( this->_print_mesh_info ) this->_mesh->print_info();

   // Print info if requested
  if( this->_print_equation_system_info ) this->_equation_system->print_info();

  return;
}

std::tr1::shared_ptr<libMesh::EquationSystems> GRINS::Simulation::get_equation_system()
{
  return _equation_system;
}

void GRINS::Simulation::check_for_restart( const GetPot& input )
{
  const std::string restart_file = input( "restart-options/restart_file", "none" );

  // Most of this was pulled from FIN-S
  if (restart_file != "none")
    {
      std::cout << " ====== Restarting from " << restart_file << std::endl;      

       // Must have correct file type to restart
      if (restart_file.rfind(".xdr") < restart_file.size())
        _equation_system->read(restart_file,libMeshEnums::DECODE,
			      //EquationSystems::READ_HEADER |  // Allow for thermochemistry upgrades
			      EquationSystems::READ_DATA |
			      EquationSystems::READ_ADDITIONAL_DATA);
      
      else if  (restart_file.rfind(".xda") < restart_file.size())
        _equation_system->read(restart_file,libMeshEnums::READ,
			      //EquationSystems::READ_HEADER |  // Allow for thermochemistry upgrades
			      EquationSystems::READ_DATA |
			      EquationSystems::READ_ADDITIONAL_DATA);

      else
        {
          std::cerr << "Error: Restart filename must have .xdr or .xda extension!" << std::endl;
          libmesh_error();
        }
      
      const std::string system_name = input("screen-options/system_name", "GRINS" );

      GRINS::MultiphysicsSystem& system = 
        _equation_system->get_system<GRINS::MultiphysicsSystem>(system_name);

      // Update the old data
      system.update();
    }

  return;
}

void GRINS::Simulation::attach_neumann_bc_funcs( std::map< std::string, GRINS::NBCContainer > neumann_bcs,
						 GRINS::MultiphysicsSystem* system )
{
  //_neumann_bc_funcs = neumann_bcs;

  if( neumann_bcs.size() > 0 )
    {
      for( std::map< std::string, GRINS::NBCContainer >::iterator bc = neumann_bcs.begin();
	   bc != neumann_bcs.end();
	   bc++ )
	{
	  std::tr1::shared_ptr<GRINS::Physics> physics = system->get_physics( bc->first );
	  physics->attach_neumann_bound_func( bc->second );
	}
    }

  return;
}

void GRINS::Simulation::attach_dirichlet_bc_funcs( std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > dbc_map,
						   GRINS::MultiphysicsSystem* system )
{
  for( std::multimap< GRINS::PhysicsName, GRINS::DBCContainer >::const_iterator it = dbc_map.begin();
       it != dbc_map.end();
       it++ )
    {
      std::tr1::shared_ptr<GRINS::Physics> physics = system->get_physics( it->first );
      
      physics->attach_dirichlet_bound_func( it->second );
    }
  return;
}

#ifdef USE_GRVY_TIMERS
void GRINS::Simulation::attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer )
{
  this->_multiphysics_system->attach_grvy_timer( grvy_timer );
  return;
}
#endif
