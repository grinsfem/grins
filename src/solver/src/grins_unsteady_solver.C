//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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


// This class
#include "grins/grins_unsteady_solver.h"

// GRINS
#include "grins/grins_enums.h"
#include "grins/solver_context.h"
#include "grins/multiphysics_sys.h"
#include "grins/time_stepping_parsing.h"
#include "grins/strategies_parsing.h"
#include "grins/solver_names.h"

// libMesh
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/getpot.h"
#include "libmesh/euler_solver.h"
#include "libmesh/euler2_solver.h"
#include "libmesh/twostep_time_solver.h"

// C++
#include <ctime>

namespace GRINS
{

  UnsteadySolver::UnsteadySolver( const GetPot& input )
    : Solver(input),
      _time_solver_name(TimeSteppingParsing::parse_time_stepper_name(input)),
      _n_timesteps( TimeSteppingParsing::parse_n_timesteps(input) ),
      _backtrack_deltat( TimeSteppingParsing::parse_backtrack_deltat(input) ),
      _theta( TimeSteppingParsing::parse_theta(input) ),
      _deltat( TimeSteppingParsing::parse_deltat(input) ),
      _target_tolerance( StrategiesParsing::parse_target_tolerance(input) ),
      _upper_tolerance( StrategiesParsing::parse_upper_tolerance(input) ),
      _max_growth( StrategiesParsing::parse_max_growth(input) ),
      _is_second_order_in_time(false)
  {
    StrategiesParsing::parse_component_norm(input,_component_norm);
  }

  void UnsteadySolver::init_time_solver(MultiphysicsSystem* system)
  {
    libMesh::UnsteadySolver* time_solver = NULL;

    if( _time_solver_name == SolverNames::libmesh_euler_solver() )
      {
        time_solver = new libMesh::EulerSolver( *(system) );

        this->set_theta<libMesh::EulerSolver>(time_solver);
      }
    else if( _time_solver_name == SolverNames::libmesh_euler2_solver() )
      {
        time_solver = new libMesh::Euler2Solver( *(system) );

        this->set_theta<libMesh::Euler2Solver>(time_solver);
      }
    else
      libmesh_error_msg("ERROR: Unsupported time stepper "+_time_solver_name);

    if (_target_tolerance)
      {
        libMesh::TwostepTimeSolver *outer_solver =
          new libMesh::TwostepTimeSolver(*system);

        outer_solver->target_tolerance = _target_tolerance;
        outer_solver->upper_tolerance = _upper_tolerance;
        outer_solver->max_growth = _max_growth;
        outer_solver->quiet = false;

        outer_solver->core_time_solver =
          libMesh::AutoPtr<libMesh::UnsteadySolver>(time_solver);
        system->time_solver = libMesh::AutoPtr<libMesh::TimeSolver>(outer_solver);
      } 
    else
      {
        system->time_solver = libMesh::AutoPtr<libMesh::TimeSolver>(time_solver);
      }

    time_solver->reduce_deltat_on_diffsolver_failure = this->_backtrack_deltat;
  }

  void UnsteadySolver::solve( SolverContext& context )
  {
    libmesh_assert( context.system );

    context.system->deltat = this->_deltat;
  
    libMesh::Real sim_time;

    if( context.output_vis ) 
      {
	context.postprocessing->update_quantities( *(context.equation_system) );
	context.vis->output( context.equation_system );
      }

    std::time_t first_wall_time = std::time(NULL);
    
    // Now we begin the timestep loop to compute the time-accurate
    // solution of the equations.
    for (unsigned int t_step=0; t_step < this->_n_timesteps; t_step++)
      {
        std::time_t latest_wall_time = std::time(NULL);

	std::cout << "==========================================================" << std::endl
		  << "   Beginning time step " << t_step  <<
                     ", t = " << context.system->time <<
                     ", dt = " << context.system->deltat <<
                     ", runtime = " << (latest_wall_time - first_wall_time) << 
                     std::endl
		  << "==========================================================" << std::endl;

        // If we have any solution-dependent Dirichlet boundaries, we
        // need to update them with the current solution.
        this->update_dirichlet_bcs(context);

	// GRVY timers contained in here (if enabled)
	context.system->solve();

	sim_time = context.system->time;

	if( context.output_vis && !((t_step+1)%context.timesteps_per_vis) )
	  {
	    context.postprocessing->update_quantities( *(context.equation_system) );
	    context.vis->output( context.equation_system, t_step, sim_time );
	  }

	if( context.output_residual && !((t_step+1)%context.timesteps_per_vis) )
	  context.vis->output_residual( context.equation_system, context.system,
                                        t_step, sim_time );

        if ( context.print_perflog && context.timesteps_per_perflog
             && !((t_step+1)%context.timesteps_per_perflog) )
          libMesh::perflog.print_log();

        if ( context.print_scalars )
          this->print_scalar_vars(context);

	// Advance to the next timestep
	context.system->time_solver->advance_timestep();
      }

    std::time_t final_wall_time = std::time(NULL);
    std::cout << "==========================================================" << std::endl
	      << "   Ending time stepping, t = " << context.system->time <<
                 ", runtime = " << (final_wall_time - first_wall_time) << 
                 std::endl
              << "==========================================================" << std::endl;


    return;
  }

  void UnsteadySolver::update_dirichlet_bcs( SolverContext& context )
  {
    // FIXME: This needs to be much more efficient and intuitive.
    // FIXME: This is only checking for nonlinear bc! This is not checking for time-dependence!
    bool have_nonlinear_dirichlet_bc = false;
    {
      const libMesh::DirichletBoundaries &db =
        *context.system->get_dof_map().get_dirichlet_boundaries();
      for (libMesh::DirichletBoundaries::const_iterator
             it = db.begin(); it != db.end(); ++it)
        {
          const libMesh::DirichletBoundary* bdy = *it;
          if (bdy->f_fem.get())
            {
              have_nonlinear_dirichlet_bc = true;
              break;
            }
        }
    }

    // Nonlinear Dirichlet constraints change as the solution does
    // FIXME: We should be updating with time-dependent BCs as well!
    if (have_nonlinear_dirichlet_bc)
      {
        context.system->reinit_constraints();
        context.system->get_dof_map().enforce_constraints_exactly(*context.system);
        context.system->get_dof_map().enforce_constraints_exactly(*context.system,
                                                                  dynamic_cast<libMesh::UnsteadySolver*>(context.system->time_solver.get())->old_local_nonlinear_solution.get());
      }
  }
} // namespace GRINS
