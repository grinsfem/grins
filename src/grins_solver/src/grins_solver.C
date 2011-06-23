//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - a low Mach number Navier-Stokes Finite-Element Solver
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
  this->_vis_output_file_prefix = input("vis-options/vis_output_file_prefix", "unknown" );
  this->_output_format          = input("vis-options/output_format", "ExodusII" );
  this->_output_vis_time_series = input("vis-options/output_vis_time_series", false);

  // Screen display options
  this->_print_mesh_info            = input("screen-options/print_mesh_info", false );
  this->_print_log_info             = input("screen-options/print_log_info", false );
  this->_solver_quiet               = input("screen-options/solver_quiet", false );
  this->_solver_verbose             = input("screen-options/solver_verbose", false );
  this->_print_equation_system_info = input("screen-options/print_equation_system_info", false );

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
  // Print mesh info if the user wants it
  if( this->_print_mesh_info ) this->_mesh->print_info();

  // Only print log info if the user requests it
  libMesh::perflog.disable_logging();
  if( this->_print_log_info ) libMesh::perflog.enable_logging();

  // Create an equation systems object.
  this->_equation_systems = new libMesh::EquationSystems(*_mesh);

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

  // Print info if requested
  if( this->_print_equation_system_info ) this->_equation_systems->print_info();

  // Get diff solver to set options
  libMesh::DiffSolver &solver = *(this->_system->time_solver->diff_solver().get());

  // Set linear/nonlinear solver options
  this->set_solver_options( solver );

  this->_system_initialized = true;

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

      // Advance to the next timestep in a transient problem
      this->_system->time_solver->advance_timestep();
    } // End time loop.

  return;
}

template< class T >
void GRINS::Solver<T>::output_visualization()
{
  this->dump_visualization( this->_vis_output_file_prefix );
  return;
}

template< class T >
void GRINS::Solver<T>::output_visualization( unsigned int time_step )
{
  std::stringstream suffix;

  suffix << time_step;

  std::string filename = this->_vis_output_file_prefix;
  filename+=suffix.str();

  this->dump_visualization( filename );

  return;
}

template< class T >
void GRINS::Solver<T>::dump_visualization( std::string filename )
{
  if( this->_vis_output_file_prefix == "unknown" )
    {
      // TODO: Need consisent way to print warning messages.
      std::cout << " WARNING in GRINS::Solver::dump_visualization :" <<
                   " using 'unknown' as file prefix since it was not set " <<
                   std::endl;
    }

  // The following is a modifed copy from the FIN-S code.
  if (this->_output_format == "tecplot" ||
      this->_output_format == "dat")
    {
      filename+=".dat";
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
  else if (this->_output_format == "tecplot_binary" ||
	   this->_output_format == "plt")
    {
      filename+=".plt";
      libMesh::TecplotIO(*(this->_mesh),true).write_equation_systems(
						filename,
						*(this->_equation_systems) );
    }
  else if (this->_output_format == "gmv")
    {
      filename+=".gmv";
      GMVIO(*(this->_mesh)).write_equation_systems(
						filename,
						*(this->_equation_systems) );
    }
  else if (this->_output_format == "vtu")
    {
      filename+=".vtu";
      VTKIO(*(this->_mesh)).write_equation_systems(
						filename,
						*(this->_equation_systems) );
    }
  else if (this->_output_format == "ExodusII")
    {
      filename+=".exo";
      ExodusII_IO(*(this->_mesh)).write_equation_systems(
						filename,
						*(this->_equation_systems) );
    }
  else if (this->_output_format.find("xda") != std::string::npos ||
	   this->_output_format.find("xdr") != std::string::npos)
    {
      filename+=this->_output_format;
      const bool binary = (this->_output_format.find("xdr") != std::string::npos);
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

  return;
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
