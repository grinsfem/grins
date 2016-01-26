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
#include "grins/unsteady_mesh_adaptive_solver.h"

// GRINS
#include "grins/common.h"
#include "grins/solver_context.h"
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/error_vector.h"

namespace GRINS
{
  UnsteadyMeshAdaptiveSolver::UnsteadyMeshAdaptiveSolver( const GetPot& input )
    : UnsteadySolver(input),
      MeshAdaptiveSolverBase( input )
  {}

  void UnsteadyMeshAdaptiveSolver::solve(  SolverContext& context )
  {
    context.system->deltat = this->_deltat;

    libMesh::Real sim_time;

    if( context.output_vis )
      {
	context.postprocessing->update_quantities( *(context.equation_system) );
	context.vis->output( context.equation_system );
      }

    // Setup MeshRefinement
    libMesh::MeshBase& mesh = context.equation_system->get_mesh();
    this->build_mesh_refinement( mesh );

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

        for ( unsigned int r_step = 0; r_step < this->_max_refinement_steps; r_step++ )
          {
            std::cout << "==========================================================" << std::endl
                      << "Adaptive Refinement Step " << r_step << std::endl
                      << "==========================================================" << std::endl;

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


             // Now we construct the data structures for the mesh refinement process
            libMesh::ErrorVector error;
            this->estimate_error_for_amr(context,error);

            // Check for convergence of error
            bool converged = this->check_for_convergence( context, error );

            if( converged )
              {
                // Break out of adaptive loop
                std::cout << "==========================================================" << std::endl
                          << "Convergence detected!" << std::endl
                          << "==========================================================" << std::endl;
                break;
              }
            else
              {
                // Only bother refining if we're not on the last step.
                if( r_step < this->_max_refinement_steps )
                  this->perform_amr(context,error);
              }

          } // End mesh adaptive loop

        // Advance to the next timestep
        context.system->time_solver->advance_timestep();

      } // End time step loop

    std::time_t final_wall_time = std::time(NULL);
    std::cout << "==========================================================" << std::endl
	      << "   Ending time stepping, t = " << context.system->time <<
                 ", runtime = " << (final_wall_time - first_wall_time) <<
                 std::endl
              << "==========================================================" << std::endl;

    // Print out the QoI, but only do it if the user asks for it
    if( context.print_qoi )
      this->print_qoi(context,std::cout);
  }

} // namespace GRINS
