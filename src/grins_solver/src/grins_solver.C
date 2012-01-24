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

#include "grins_solver.h"

// For instantiating systems
#include "multiphysics_sys.h"

// libMesh I/O classes
#include "gmv_io.h"
#include "tecplot_io.h"
#include "exodusII_io.h"
#include "vtk_io.h"

#include <iostream>

template< class T >
GRINS::Solver<T>::Solver( )
  : _system_initialized(false)
{
  return;
}

template< class T >
GRINS::Solver<T>::~Solver()
{
  // If we initialized the system, we need to destroy the
  // EquationSystems object.
  if(this->_system_initialized)
    {
      delete this->_equation_systems;
    }
  return;
}

template< class T >
void GRINS::Solver<T>::read_input_options( const GetPot& input )
{
  // Linear/Nonlinear solver options
  this->_max_nonlinear_iterations    = input("linear-nonlinear-solver/max_nonlinear_iterations", 10 );
  this->_relative_step_tolerance     = input("linear-nonlinear-solver/relative_step_tolerance", 1.e-6 );
  this->_absolute_step_tolerance     = input("linear-nonlinear-solver/absolute_step_tolerance", 0.0 );
  this->_relative_residual_tolerance = input("linear-nonlinear-solver/relative_residual_tolerance", 1.e-16 );
  this->_absolute_residual_tolerance = input("linear-nonlinear-solver/absolute_residual_tolerance", 0.0 );
  this->_max_linear_iterations       = input("linear-nonlinear-solver/max_linear_iterations", 500 );
  this->_initial_linear_tolerance    = input("linear-nonlinear-solver/initial_linear_tolerance", 1.e-3 );

  // Unsteady solver options
  this->_transient   = input("unsteady-solver/transient", false );
  this->_theta       = input("unsteady-solver/theta", 0.5 );
  this->_n_timesteps = input("unsteady-solver/n_timesteps", 1 );
  this->_deltat      = input("unsteady-solver/deltat", 0.0 ); //TODO: Better default here?

  // Visualization options
  this->_vis_output_file_prefix   = input("vis-options/vis_output_file_prefix", "unknown" );
  this->_output_vis_time_series   = input("vis-options/output_vis_time_series", false);
  this->_output_residual          = input("vis-options/output_residual", false );
  this->_output_unsteady_residual = input("vis-options/output_unsteady_residual", false );

  unsigned int num_formats = input.vector_variable_size("vis-options/output_format");

  // If no format specified, default to ExodusII only
  if( num_formats == 0 )
    {
      _output_format.push_back("ExodusII");
    }

  for( unsigned int i = 0; i < num_formats; i++ )
    {
      _output_format.push_back( input("vis-options/output_format", "DIE", i ) );
    }
  
  // Screen display options
  this->_solver_quiet               = input("screen-options/solver_quiet", false );
  this->_solver_verbose             = input("screen-options/solver_verbose", false );
  

  return;
}

template< class T >
libMesh::Mesh* GRINS::Solver<T>::get_mesh()
{
  return this->_mesh;
}

template< class T >
void GRINS::Solver<T>::set_mesh( libMesh::Mesh *mesh )
{
  this->_mesh = mesh;
  return;
}

template< class T >
void GRINS::Solver<T>::set_solver_options( libMesh::DiffSolver& solver  )
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

template< class T >
void GRINS::Solver<T>::initialize_system( std::string system_name, GetPot& input )
{
  // Declare the system and its variables.
  libMesh::EquationSystems *es = this->_equation_systems;
  this->_system = &es->add_system<T> (system_name);

  this->_system->read_input_options( input );

  // Solve this as a time-dependent or steady system
  if (this->_transient)
    {
      // Setup time solver
      libMesh::EulerSolver* time_solver = new libMesh::EulerSolver( *(this->_system) );

      this->_system->time_solver = AutoPtr<TimeSolver>(time_solver);

      // Set theta parameter for time-stepping scheme
      time_solver->theta = this->_theta;
    }
  else
    {
      libMesh::SteadySolver* time_solver = new libMesh::SteadySolver( *(this->_system) );
      this->_system->time_solver = AutoPtr<TimeSolver>(time_solver);
      libmesh_assert( this->_n_timesteps == 1 );
    }

  // Initialize the system
  this->_equation_systems->init();

  // Get diff solver to set options
  libMesh::DiffSolver &solver = *(this->_system->time_solver->diff_solver().get());

  // Set linear/nonlinear solver options
  this->set_solver_options( solver );

  return;
}

template< class T >
void GRINS::Solver<T>::solve()
{
  this->_system->deltat = this->_deltat;

  // Now we begin the timestep loop to compute the time-accurate
  // solution of the equations.
  for (unsigned int t_step=0; t_step < this->_n_timesteps; t_step++)
    {
      // GRVY timers contained in here (if enabled)
      this->_system->solve();

      // Dump out time series visualization if user wants it.
      if( this->_output_vis_time_series )
	{
	  this->output_visualization( t_step );
	}

      if( this->_output_residual )
	{
	  this->output_residual_vis( t_step );
	}

      if( this->_output_unsteady_residual )
	{
	  this->output_unsteady_residual_vis( t_step );
	}

      // Advance to the next timestep in a transient problem
      this->_system->time_solver->advance_timestep();
    } // End time loop.

  return;
}

template< class T >
void GRINS::Solver<T>::output_visualization()
{
  this->dump_visualization( this->_vis_output_file_prefix, 0 );

  return;
}

template< class T >
void GRINS::Solver<T>::output_visualization( unsigned int time_step )
{
  std::stringstream suffix;

  suffix << time_step;

  std::string filename = this->_vis_output_file_prefix;
  filename+="."+suffix.str();

  this->dump_visualization( filename, time_step );

  return;
}

template< class T >
void GRINS::Solver<T>::output_residual_vis( const unsigned int time_step )
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

template< class T >
void GRINS::Solver<T>::output_unsteady_residual_vis( const unsigned int time_step )
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

template< class T >
void GRINS::Solver<T>::dump_visualization( const std::string filename_prefix, const int time_step )
{
  if( this->_vis_output_file_prefix == "unknown" )
    {
      // TODO: Need consisent way to print warning messages.
      std::cout << " WARNING in GRINS::Solver::dump_visualization :" <<
	" using 'unknown' as file prefix since it was not set " <<
	std::endl;
    }
  
  for( std::vector<std::string>::const_iterator format = _output_format.begin();
       format != _output_format.end();
       format ++ )
    {
      // The following is a modifed copy from the FIN-S code.
      if ((*format) == "tecplot" ||
	  (*format) == "dat")
	{
	  std::string filename = filename_prefix+".dat";
	  libMesh::TecplotIO(*(this->_mesh),false).write_equation_systems(
									  filename,
									  *(this->_equation_systems) );
	  
	  // Left this here as an example from FIN-S if we need to
	  // handle boundary meshes separately in the future.
	  /*
	    if (have_boundary_data)
	    {
	    sprintf( filechar, "%s-surf-%05d.dat",
	    output_name.c_str(),
	    write_soln_number );
	    
	    libMesh::TecplotIO(*boundary_mesh,false).write_equation_systems(
	    std::string(filechar),
	    *_boundary_equation_systems );
	    }
	  */
	}
      else if ((*format) == "tecplot_binary" ||
	       (*format) == "plt")
	{
	  std::string filename = filename_prefix+".plt";
	  libMesh::TecplotIO(*(this->_mesh),true).write_equation_systems(
									 filename,
									 *(this->_equation_systems) );
	}
      else if ((*format) == "gmv")
	{
	  std::string filename = filename_prefix+".gmv";
	  GMVIO(*(this->_mesh)).write_equation_systems(
						       filename,
						       *(this->_equation_systems) );
	}
      else if ((*format) == "vtu")
	{
	  std::string filename = filename_prefix+".vtu";
	  VTKIO(*(this->_mesh)).write_equation_systems(
						       filename,
						       *(this->_equation_systems) );
	}
      else if ((*format) == "ExodusII")
	{
	  std::string filename = filename_prefix+".exo";
	  
	  // The "1" is hardcoded for the number of time steps because the ExodusII manual states that
	  // it should be the number of timesteps within the file. Here, we are explicitly only doing 
	  // one timestep per file.
	  ExodusII_IO(*(this->_mesh)).write_timestep(
						     filename,
						     *(this->_equation_systems),
						     1,
						     this->_system->time );
	}
      else if ((*format).find("xda") != std::string::npos ||
	       (*format).find("xdr") != std::string::npos)
	{
	  std::string filename = filename_prefix+"."+(*format);
	  const bool binary = ((*format).find("xdr") != std::string::npos);
	  (this->_equation_systems)->write( filename,
					    binary ? libMeshEnums::ENCODE : libMeshEnums::WRITE,
					    EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA );
	}
      else
	{
	  // TODO: Do we want to use this to error throughout the code?
	  // TODO: (at least need to pass/print some message/string) - sahni
	  libmesh_error();
	}
    } // End loop over formats

  return;
}

template< class T >
T* GRINS::Solver<T>::get_system()
{
  return this->_system;
}


#ifdef USE_GRVY_TIMERS
template< class T >
void GRINS::Solver<T>::attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer )
{
  _timer = grvy_timer;

  // Attach timer to system
  this->_system->attach_grvy_timer( grvy_timer );

  return;
}
#endif

// Instantiate class
template class GRINS::Solver<GRINS::MultiphysicsSystem>;
