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

// libMesh
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/getpot.h"
#include "libmesh/euler_solver.h"
#include "libmesh/twostep_time_solver.h"

// C++
#include <ctime>

namespace GRINS
{

  UnsteadySolver::UnsteadySolver( const GetPot& input )
    : Solver(input),
      _n_timesteps( input("unsteady-solver/n_timesteps", 1 ) ),
      _backtrack_deltat( input("unsteady-solver/backtrack_deltat", 0 ) ),
      _theta( input("unsteady-solver/theta", 0.5 ) ),
      /*! \todo Is this the best default for delta t?*/
      _deltat( input("unsteady-solver/deltat", 0.0 ) ),
      _target_tolerance( input("unsteady-solver/target_tolerance", 0.0 ) ),
      _upper_tolerance( input("unsteady-solver/upper_tolerance", 0.0 ) ),
      _max_growth( input("unsteady-solver/max_growth", 0.0 ) )
  {
    const unsigned int n_component_norm =
      input.vector_variable_size("unsteady-solver/component_norm");
    for (unsigned int i=0; i != n_component_norm; ++i)
      {
        const std::string current_norm = input("component_norm", std::string("L2"), i);
        // TODO: replace this with string_to_enum with newer libMesh
        if (current_norm == "GRINSEnums::L2")
          _component_norm.set_type(i, libMesh::L2);
        else if (current_norm == "GRINSEnums::H1")
          _component_norm.set_type(i, libMesh::H1);
        else
          libmesh_not_implemented();
      }


  }

  UnsteadySolver::~UnsteadySolver()
  {
    return;
  }

  void UnsteadySolver::init_time_solver(MultiphysicsSystem* system)
  {
    libMesh::EulerSolver* time_solver = new libMesh::EulerSolver( *(system) );

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

    // Set theta parameter for time-stepping scheme
    time_solver->theta = this->_theta;
    time_solver->reduce_deltat_on_diffsolver_failure = this->_backtrack_deltat;

    return;
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
        //
        // FIXME: This needs to be much more efficient and intuitive.
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
        if (have_nonlinear_dirichlet_bc)
          {
            context.system->reinit_constraints();
            context.system->get_dof_map().enforce_constraints_exactly(*context.system);
            context.system->get_dof_map().enforce_constraints_exactly(*context.system,
                                                                      dynamic_cast<libMesh::UnsteadySolver*>(context.system->time_solver.get())->old_local_nonlinear_solution.get());
          }

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
          for (unsigned int v=0; v != context.system->n_vars(); ++v)
            if (context.system->variable(v).type().family ==
                libMesh::SCALAR)
              {
                std::cout << context.system->variable_name(v) <<
                             " = {";
                std::vector<libMesh::dof_id_type> scalar_indices;
                context.system->get_dof_map().SCALAR_dof_indices
                  (scalar_indices, v);
                if (scalar_indices.size())
                  std::cout <<
                    context.system->current_solution(scalar_indices[0]);
                for (unsigned int i=1; i < scalar_indices.size();
                     ++i)
                  std::cout << ", " <<
                    context.system->current_solution(scalar_indices[i]);
                std::cout << '}' << std::endl;
              }

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

} // namespace GRINS
