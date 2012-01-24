//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "simulation.h"

Simulation::Simulation( const GetPot& input,
			GRINS::PhysicsFactory* physics_factory,
			GRINS::MeshBuilder* mesh_builder,
			GRINS::SolverFactory* solver_factory,
			GRINS::VisualizationFactory* vis_factory )
  :  _mesh( mesh_builder->build() ),
     _equation_system( new libMesh::EquationSystems( *_mesh ) ),
     _solver( solver_factory->build() ),
     _vis( vis_factory->build() )
{
  // Read input options relevant for this object
  this->read_input_options( input );
  
  // Only print libMesh logging info if the user requests it
  libMesh::perflog.disable_logging();
  if( this->_print_log_info ) libMesh::perflog.enable_logging();

  _solver->initialize( input, _equation_system, physics_factory->build() );

  return;
}

Simulation::~Simulation()
{
  return;
}

void Simulation::run()
{
  this->print_sim_info();
  
  _solver->solve();

  return;
}

void Simulation::output_vis()
{
  _vis->output( _equation_system );
  return;
}

void Simulation::read_input_options( const GetPot& input )
{
  this->_print_mesh_info            = input("screen-options/print_mesh_info", false );
  this->_print_log_info             = input("screen-options/print_log_info", false );
  this->_print_equation_system_info = input("screen-options/print_equation_system_info", false );

  return;
}

void Simulation::print_sim_info()
{
  // Print mesh info if the user wants it
  if( this->_print_mesh_info ) this->_mesh->print_info();

   // Print info if requested
  if( this->_print_equation_system_info ) this->_equation_systems->print_info();

  return;
}
