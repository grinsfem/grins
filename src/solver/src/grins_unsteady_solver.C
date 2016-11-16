//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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
#include "libmesh/function_base.h"
#include "libmesh/twostep_time_solver.h"
#include "libmesh/newmark_solver.h"
#include "libmesh/function_base.h"

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
      _adapt_time_step_options(input),
      _is_second_order_in_time(false)
  {}

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
    else if( _time_solver_name == SolverNames::libmesh_newmark_solver() )
      {
        time_solver = new libMesh::NewmarkSolver( *(system) );
        _is_second_order_in_time = true;
      }
    else
      libmesh_error_msg("ERROR: Unsupported time stepper "+_time_solver_name);

    if( _adapt_time_step_options.is_time_adaptive() )
      {
        libMesh::TwostepTimeSolver *outer_solver =
          new libMesh::TwostepTimeSolver(*system);

        outer_solver->target_tolerance = _adapt_time_step_options.target_tolerance();
        outer_solver->upper_tolerance = _adapt_time_step_options.upper_tolerance();
        outer_solver->max_growth = _adapt_time_step_options.max_growth();
        outer_solver->component_norm = _adapt_time_step_options.component_norm();
        outer_solver->quiet = false;

        outer_solver->core_time_solver =
          libMesh::UniquePtr<libMesh::UnsteadySolver>(time_solver);
        system->time_solver = libMesh::UniquePtr<libMesh::TimeSolver>(outer_solver);
      }
    else
      {
        system->time_solver = libMesh::UniquePtr<libMesh::TimeSolver>(time_solver);
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

    // We may need to initialize acceleration for second order solvers
    if( _is_second_order_in_time )
      this->init_second_order_in_time_solvers(context);

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
    bool have_nonlinear_dirichlet_bc = false;
    bool have_time_dependence = false;
    {
      const libMesh::DirichletBoundaries &db =
        *context.system->get_dof_map().get_dirichlet_boundaries();

      for (libMesh::DirichletBoundaries::const_iterator
             it = db.begin(); it != db.end(); ++it)
        {
          const libMesh::DirichletBoundary* bdy = *it;

          // If we have a FEMFunctionBase, we assume nonlinearity
          if (bdy->f_fem.get())
              have_nonlinear_dirichlet_bc = true;

          // Check for time-dependence of FunctionBase
          if( bdy->f.get() )
              if( bdy->f->is_time_dependent() )
                  have_time_dependence = true;

          if( have_nonlinear_dirichlet_bc || have_time_dependence )
            break;

        } // End loop over DirichletBoundaries
    }


    // Nonlinear Dirichlet constraints change as the solution does
    // and time-dependent constraints have to be updated
    if (have_nonlinear_dirichlet_bc || have_time_dependence )
      {
        context.system->reinit_constraints();
        context.system->get_dof_map().enforce_constraints_exactly(*context.system);
        context.system->get_dof_map().enforce_constraints_exactly(*context.system,
                                                                  dynamic_cast<libMesh::UnsteadySolver*>(context.system->time_solver.get())->old_local_nonlinear_solution.get());
      }
  }

  void UnsteadySolver::init_second_order_in_time_solvers( SolverContext& context )
  {
    // Right now, only Newmark is available so we cast directly to that
    libMesh::TimeSolver& base_time_solver = context.system->get_time_solver();

    libMesh::NewmarkSolver& time_solver = libMesh::cast_ref<libMesh::NewmarkSolver&>(base_time_solver);

    // If there's a restart, the acceleration should already be there
    if( context.have_restart )
      time_solver.set_initial_accel_avail(true);

    // Otherwise, we need to compute it
    else
      {
        libMesh::out << "==========================================================" << std::endl
                     << "            Computing Initital Acceleration" << std::endl
                     << "==========================================================" << std::endl;

        time_solver.compute_initial_accel();
      }
  }

} // namespace GRINS
