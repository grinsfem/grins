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
#include "grins/steady_mesh_adaptive_solver.h"

// GRINS
#include "grins/solver_context.h"
#include "grins/multiphysics_sys.h"
#include "grins/composite_qoi.h"
#include "grins/common.h"

// libMesh
#include "libmesh/error_vector.h"
#include "libmesh/steady_solver.h"
#include "libmesh/adjoint_refinement_estimator.h"

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

    /*! \todo This output cannot be toggled in the input file, but it should be able to be. */
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

        if( context.output_vis )
          {
            context.postprocessing->update_quantities( *(context.equation_system) );
            context.vis->output( context.equation_system );
          }

        // Solve adjoint system
        if(context.do_adjoint_solve)
          this->steady_adjoint_solve(context);

        if(context.output_adjoint)
          context.vis->output_adjoint(context.equation_system, context.system);

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

	// Get the global error estimate if you can and are asked to
	if( this->compute_QoI_error_estimate )
	{
	  // First check if we are using the Adjoint Refinement Error Estimator
	  if(context.error_estimator->type() == libMesh::ADJOINT_REFINEMENT)
	  {
	    // Loop over QoIs and print error error estimates
	    for(unsigned int i = 0; i != context.system->qoi.size(); i++)
	    {
	      libMesh::AdjointRefinementEstimator* adjoint_ref_error_estimator = libMesh::libmesh_cast_ptr<libMesh::AdjointRefinementEstimator*>( context.error_estimator.get() );
	      std::cout<<"The error estimate for QoI("<<i<<") is: "<<adjoint_ref_error_estimator->get_global_QoI_error_estimate(i)<<std::endl;
	    }

	  }
	  else
	  {
	    std::string warning = "You asked for QoI error estimates but did not use an Adjoint Refinement Error Estimator!\n";
	    warning += "Please use the ADJOINT_REFINEMENT option for the estimator_type if you want QoI error estimates.\n";
	    grins_warning(warning);
	  }
	}

        // Check for convergence of error
        std::cout << "==========================================================" << std::endl
                  << "Checking convergence" << std::endl
                  << "==========================================================" << std::endl;
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

                // It's helpful to print the qoi along the way, but only do it if the user
                // asks for it
                if( context.print_qoi )
                  {
                    context.system->assemble_qoi();
                    const CompositeQoI* my_qoi = libMesh::libmesh_cast_ptr<const CompositeQoI*>(context.system->get_qoi());
                    my_qoi->output_qoi( std::cout );
                    std::cout << std::endl;
                  }
              }
          }

      } // r_step for-loop

    return;
  }

  void SteadyMeshAdaptiveSolver::adjoint_qoi_parameter_sensitivity
    (SolverContext& context,
     const libMesh::QoISet&          qoi_indices,
     const libMesh::ParameterVector& parameters_in,
     libMesh::SensitivityData&       sensitivities) const
  {
    context.system->adjoint_qoi_parameter_sensitivity
      (qoi_indices, parameters_in, sensitivities);
  }

  void SteadyMeshAdaptiveSolver::forward_qoi_parameter_sensitivity
    (SolverContext& context,
     const libMesh::QoISet&          qoi_indices,
     const libMesh::ParameterVector& parameters_in,
     libMesh::SensitivityData&       sensitivities) const
  {
    context.system->forward_qoi_parameter_sensitivity
      (qoi_indices, parameters_in, sensitivities);

    if( context.output_residual_sensitivities )
      context.vis->output_residual_sensitivities
        ( context.equation_system, context.system, parameters_in );

    if( context.output_solution_sensitivities )
      context.vis->output_solution_sensitivities
        ( context.equation_system, context.system, parameters_in );
  }


} // end namespace GRINS
