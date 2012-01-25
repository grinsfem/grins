//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010,2011 The PECOS Development Team
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

#include <iostream>

#include "grins_solver.h"

GRINS::Solver::Solver( const GetPot& input )
  : _system(NULL),
    _max_nonlinear_iterations( input("linear-nonlinear-solver/max_nonlinear_iterations", 10 ) ),
    _relative_step_tolerance( input("linear-nonlinear-solver/relative_step_tolerance", 1.e-6 ) ),
    _absolute_step_tolerance( input("linear-nonlinear-solver/absolute_step_tolerance", 0.0 ) ),
    _relative_residual_tolerance( input("linear-nonlinear-solver/relative_residual_tolerance", 1.e-16 ) ),
    _absolute_residual_tolerance( input("linear-nonlinear-solver/absolute_residual_tolerance", 0.0 ) ),
    _max_linear_iterations( input("linear-nonlinear-solver/max_linear_iterations", 500 ) ),
    _initial_linear_tolerance( input("linear-nonlinear-solver/initial_linear_tolerance", 1.e-3 ) ),
    _solver_quiet( input("screen-options/solver_quiet", false ) ),
    _solver_verbose( input("screen-options/solver_verbose", false ) ),
    _system_name( input("screen-options/system_name", "GRINS" ) )
{
  return;
}


GRINS::Solver::~Solver()
{
  return;
}

void GRINS::Solver::initialize( GetPot& input, 
				libMesh::EquationSystems equation_system,
				GRINS::PhysicsList& physics_list )
{
  // Declare the system and its variables.
  _system = equation_system->add_system<GRINS::MultiphysicsSystem>( _system_name );

  _system->attach_physics_list( physics_list );

  _system->read_input_options( input );

  // Defined in subclasses depending on the solver used.
  this->init_time_solver();

  // Initialize the system
  this->_equation_systems->init();

  // Get diff solver to set options
  libMesh::DiffSolver &solver = *(this->_system->time_solver->diff_solver().get());

  // Set linear/nonlinear solver options
  this->set_solver_options( solver );

  return;
}

void GRINS::Solver::output_residual_vis( const unsigned int time_step )
{
  std::stringstream suffix;
  suffix << time_step;

  std::string filename = this->_vis_output_file_prefix+"_residual";

  filename+="."+suffix.str();

  // Idea is that this->rhs stashes the residual. Thus, when we swap
  // with the solution, we should be dumping the residual. Then, we swap
  // back once we're done outputting.

  // Swap solution with computed residual
  this->_system->solution->swap( *(this->_system->rhs) );
  this->_equation_systems->update();
  
  this->dump_visualization( filename, time_step );
  
  // Now swap back and reupdate
  this->_system->solution->swap( *(this->_system->rhs) );
  this->_equation_systems->update();

  return;
}


void GRINS::Solver::output_unsteady_residual_vis( const unsigned int time_step )
{
  if( !this->_transient  &&
      this->_output_unsteady_residual )
    {
      std::cerr << "WARNING: Doesn't make sense to output unsteady residual" 
		<< " for a steady problem. Not producing output."
		<< std::endl;
      return;
    }

  std::stringstream suffix;
  suffix << time_step;

  std::string filename = this->_vis_output_file_prefix+"_unsteady_residual";

  filename+="."+suffix.str();

  // For the unsteady residual, we just want to evaluate F(u) from
  // dU/dt = F(u). What we do is swap out the time solver to a
  // SteadySolver and reassemble the residual. Then, we'll need to swap
  // the solution and the rhs vector stashed in the system. Once we're done,
  // we'll reset the time solver pointer back to the original guy.
  
  AutoPtr<TimeSolver> prev_time_solver(this->_system->time_solver);

  libMesh::SteadySolver* steady_solver = new libMesh::SteadySolver( *(this->_system) );

  this->_system->time_solver = AutoPtr<TimeSolver>(steady_solver);

  this->_system->assembly( true /*residual*/, false /*jacobian*/ );
  this->_system->rhs->close();

  // Swap solution with newly computed residual
  this->_system->solution->swap( *(this->_system->rhs) );
  // Update equation systems
  this->_equation_systems->update();
  
  this->dump_visualization( filename, time_step );
  
  // Now swap back and reupdate
  this->_system->solution->swap( *(this->_system->rhs) );
  this->_equation_systems->update();

  this->_system->time_solver = prev_time_solver;

  return;
}

void GRINS::Solver::set_solver_options( libMesh::DiffSolver& solver  )
{
  solver.quiet                       = this->_solver_quiet;
  solver.verbose                     = this->_solver_verbose;
  solver.max_nonlinear_iterations    = this->_max_nonlinear_iterations;
  solver.relative_step_tolerance     = this->_relative_step_tolerance;
  solver.absolute_step_tolerance     = this->_absolute_step_tolerance;
  solver.relative_residual_tolerance = this->_relative_residual_tolerance;
  solver.absolute_residual_tolerance = this->_absolute_residual_tolerance;
  solver.max_linear_iterations       = this->_max_linear_iterations;
  solver.initial_linear_tolerance    = this->_initial_linear_tolerance;

  return;
}

#ifdef USE_GRVY_TIMERS
void GRINS::Solver::attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer )
{
  _timer = grvy_timer;

  // Attach timer to system
  this->_system->attach_grvy_timer( grvy_timer );

  return;
}
#endif

// Instantiate class
template class GRINS::Solver<GRINS::MultiphysicsSystem>;
