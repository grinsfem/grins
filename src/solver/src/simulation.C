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

// This class
#include "grins/simulation.h"

// GRINS
#include "grins/simulation_builder.h"
#include "grins/multiphysics_sys.h"
#include "grins/solver_context.h"

namespace GRINS
{

  Simulation::Simulation( const GetPot& input,
			  SimulationBuilder& sim_builder )
    :  _mesh( sim_builder.build_mesh(input) ),
       _equation_system( new libMesh::EquationSystems( *_mesh ) ),
       _solver( sim_builder.build_solver(input) ),
       _system_name( input("screen-options/system_name", "GRINS" ) ),
       _multiphysics_system( &(_equation_system->add_system<MultiphysicsSystem>( _system_name )) ),
       _vis( sim_builder.build_vis(input) ),
       _qoi( sim_builder.build_qoi(input) ),
       _postprocessing( sim_builder.build_postprocessing(input) ),
       _print_mesh_info( input("screen-options/print_mesh_info", false ) ),
       _print_log_info( input("screen-options/print_log_info", false ) ),
       _print_equation_system_info( input("screen-options/print_equation_system_info", false ) ),
       _print_qoi( input("screen-options/print_qoi", false ) ),
       _output_vis( input("vis-options/output_vis", false ) ),
       _output_residual( input( "vis-options/output_residual", false ) )
  {
    // Only print libMesh logging info if the user requests it
    libMesh::perflog.disable_logging();
    if( this->_print_log_info ) libMesh::perflog.enable_logging();

    PhysicsList physics_list = sim_builder.build_physics(input);

    _multiphysics_system->attach_physics_list( physics_list );

    _multiphysics_system->read_input_options( input );

    // This *must* be done before equation_system->init
    this->attach_dirichlet_bc_funcs( sim_builder.build_dirichlet_bcs(), _multiphysics_system );

    /* Postprocessing needs to be initialized before the solver since that's 
       where equation_system gets init'ed */
    _postprocessing->initialize( *_multiphysics_system, *_equation_system );

    _solver->initialize( input, _equation_system, _multiphysics_system );

    // This *must* be done after equation_system->init in order to get variable indices
    this->attach_neumann_bc_funcs( sim_builder.build_neumann_bcs( *_equation_system ), _multiphysics_system );

    // If the user actually asks for a QoI, then we add it.
    if( this->_qoi.use_count() > 0 )
      {
	// This *must* be done after equation_system->init in order to get variable indices
	this->_qoi->init(input, *_multiphysics_system );
      
	/*! \todo We're missing the qoi's init_context call by putting it after equation_system->init,
	  but we also need to be able to get system variable numbers... */
	/* Note that we are effectively transfering ownership of the qoi pointer because
	   it will be cloned in _multiphysics_system and all the calculations are done there. */
	
	_multiphysics_system->attach_qoi( &(*(this->_qoi)) );
      }

    this->check_for_restart( input );

    return;
  }

  Simulation::~Simulation()
  {
    return;
  }

  void Simulation::run()
  {
    this->print_sim_info();

    SolverContext context;
    context.system = _multiphysics_system;
    context.equation_system = _equation_system;
    context.vis = _vis;
    context.output_vis = _output_vis;
    context.output_residual = _output_residual;
    context.postprocessing = _postprocessing;

    _solver->solve( context );

    if( this->_print_qoi )
      {
	_multiphysics_system->assemble_qoi( libMesh::QoISet( *_multiphysics_system ) );
	const QoIBase* my_qoi = libmesh_cast_ptr<const QoIBase*>(this->_multiphysics_system->get_qoi());
	my_qoi->output_qoi( std::cout );
      }

    return;
  }

  void Simulation::print_sim_info()
  {
    // Print mesh info if the user wants it
    if( this->_print_mesh_info ) this->_mesh->print_info();

    // Print info if requested
    if( this->_print_equation_system_info ) this->_equation_system->print_info();

    return;
  }

  std::tr1::shared_ptr<libMesh::EquationSystems> Simulation::get_equation_system()
  {
    return _equation_system;
  }

  Number Simulation::get_qoi( unsigned int qoi_index ) const
  {
    const QoIBase* qoi = libmesh_cast_ptr<const QoIBase*>(this->_multiphysics_system->get_qoi());
    return qoi->get_qoi(qoi_index);
  }

  void Simulation::check_for_restart( const GetPot& input )
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

	MultiphysicsSystem& system = 
	  _equation_system->get_system<MultiphysicsSystem>(system_name);

	// Update the old data
	system.update();
      }

    return;
  }

  void Simulation::attach_neumann_bc_funcs( std::map< std::string, NBCContainer > neumann_bcs,
					    MultiphysicsSystem* system )
  {
    //_neumann_bc_funcs = neumann_bcs;

    if( neumann_bcs.size() > 0 )
      {
	for( std::map< std::string, NBCContainer >::iterator bc = neumann_bcs.begin();
	     bc != neumann_bcs.end();
	     bc++ )
	  {
	    std::tr1::shared_ptr<Physics> physics = system->get_physics( bc->first );
	    physics->attach_neumann_bound_func( bc->second );
	  }
      }

    return;
  }

  void Simulation::attach_dirichlet_bc_funcs( std::multimap< PhysicsName, DBCContainer > dbc_map,
					      MultiphysicsSystem* system )
  {
    for( std::multimap< PhysicsName, DBCContainer >::const_iterator it = dbc_map.begin();
	 it != dbc_map.end();
	 it++ )
      {
	std::tr1::shared_ptr<Physics> physics = system->get_physics( it->first );
      
	physics->attach_dirichlet_bound_func( it->second );
      }
    return;
  }

#ifdef GRINS_USE_GRVY_TIMERS
  void Simulation::attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer )
  {
    this->_multiphysics_system->attach_grvy_timer( grvy_timer );
    return;
  }
#endif

} // namespace GRINS
