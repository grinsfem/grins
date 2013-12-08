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

// This class
#include "grins/steady_mesh_adaptive_solver.h"

// GRINS
#include "grins/solver_context.h"
#include "grins/multiphysics_sys.h"
#include "grins/composite_qoi.h"

// libMesh
#include "libmesh/error_vector.h"
#include "libmesh/steady_solver.h"

namespace GRINS
{

  SteadyMeshAdaptiveSolver::SteadyMeshAdaptiveSolver( const GetPot& input )
    : MeshAdaptiveSolverBase( input )
  {
    return;
  }

  SteadyMeshAdaptiveSolver::~SteadyMeshAdaptiveSolver()
  {
    return;
  }

  void SteadyMeshAdaptiveSolver::init_time_solver( MultiphysicsSystem* system )
  {
    libMesh::SteadySolver* time_solver = new libMesh::SteadySolver( *(system) );

    system->time_solver = libMesh::AutoPtr<libMesh::TimeSolver>( time_solver );

    return;
  }

  void SteadyMeshAdaptiveSolver::solve( SolverContext& context )
  {
    // Mesh and mesh refinement
    libMesh::MeshBase& mesh = context.equation_system->get_mesh();
    this->build_mesh_refinement( mesh );

    // This output cannot be toggled in the input file.
    std::cout << "==========================================================" << std::endl
              << "Performing " << this->_max_refinement_steps << " adaptive refinements" << std::endl
              << "==========================================================" << std::endl;

    // GRVY timers contained in here (if enabled)
    for ( unsigned int r_step = 0; r_step < this->_max_refinement_steps; r_step++ )
      {
        std::cout << "==========================================================" << std::endl
                  << "Adaptive Refinement Step " << r_step << std::endl
                  << "==========================================================" << std::endl;

        // Solve the forward problem
        context.system->solve();

        libMesh::NumericVector<Number>& primal_solution = *(context.system->solution);
        if( context.output_vis )
          {
            context.postprocessing->update_quantities( *(context.equation_system) );
            context.vis->output( context.equation_system );
          }

        // Solve adjoint system
        if( _do_adjoint_solve )
          {
            std::cout << "==========================================================" << std::endl
                      << "Solving adjoint problem." << std::endl
                      << "==========================================================" << std::endl;
            context.system->adjoint_solve();
            context.system->set_adjoint_already_solved(true);
          }

        // At the moment output data is overwritten every mesh refinement step
        if( context.output_vis && this->_output_adjoint_sol && _do_adjoint_solve )
          {
            libMesh::NumericVector<Number>& dual_solution = context.system->get_adjoint_solution(0);

            // Swap primal and dual to write out dual solution
            primal_solution.swap( dual_solution );          
            context.vis->output( context.equation_system );
            primal_solution.swap( dual_solution );          
          }

        if( context.output_residual )
          {
            context.vis->output_residual( context.equation_system, context.system );
          }

        // Now we construct the data structures for the mesh refinement process 
        libMesh::ErrorVector error;
        
        std::cout << "==========================================================" << std::endl
                  << "Estimating error" << std::endl
                  << "==========================================================" << std::endl;
        context.error_estimator->estimate_error( *context.system, error );

        // Plot error vector
        if( this->_plot_cell_errors )
          {
            error.plot_error( this->_error_plot_prefix+".exo", mesh );
          }

        // Check for convergence of error
        std::cout << "==========================================================" << std::endl
                  << "Checking convergence" << std::endl
                  << "==========================================================" << std::endl;
        bool converged = this->check_for_convergence( error );
        
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
            // Only bother refining if we're on the last step.
            if( r_step < this->_max_refinement_steps -1 )
              {
                std::cout << "==========================================================" << std::endl
                          << "Performing Mesh Refinement" << std::endl
                          << "==========================================================" << std::endl;

                this->flag_elements_for_refinement( error );
                _mesh_refinement->refine_and_coarsen_elements();
    
                // Dont forget to reinit the system after each adaptive refinement!
                context.equation_system->reinit();

                // This output cannot be toggled in the input file.
                std::cout << "==========================================================" << std::endl
                          << "Refined mesh to " << std::setw(12) << mesh.n_active_elem() 
                          << " active elements" << std::endl
                          << "            " << std::setw(16) << context.system->n_active_dofs() 
                          << " active dofs" << std::endl
                          << "==========================================================" << std::endl;

                context.system->assemble_qoi();
                const CompositeQoI* my_qoi = libmesh_cast_ptr<const CompositeQoI*>(context.system->get_qoi());
                my_qoi->output_qoi( std::cout );
                std::cout << std::endl;
              }
          }

      } // r_step for-loop

    return;
  }

} // end namespace GRINS
