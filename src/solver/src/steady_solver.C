//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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
#include "grins/steady_solver.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/solver_context.h"

// libMesh
#include "libmesh/auto_ptr.h"
#include "libmesh/dof_map.h"
#include "libmesh/getpot.h"
#include "libmesh/steady_solver.h"
#include "libmesh/linear_solver.h"

namespace GRINS
{

  SteadySolver::SteadySolver( const GetPot& input )
    : Solver( input )
  {
    return;
  }

  SteadySolver::~SteadySolver()
  {
    return;
  }

  void SteadySolver::init_time_solver(MultiphysicsSystem* system)
  {
    libMesh::SteadySolver* time_solver = new libMesh::SteadySolver( *(system) );

    system->time_solver = std::unique_ptr<libMesh::TimeSolver>(time_solver);
    return;
  }

  void SteadySolver::solve( SolverContext& context )
  {
    libmesh_assert( context.system );

    if( context.output_vis )
      {
        context.postprocessing->update_quantities( *(context.equation_system) );
        context.vis->output( context.equation_system );
      }

    context.system->solve();

    if ( context.print_scalars )
      this->print_scalar_vars(context);

    if( context.do_adjoint_solve )
      {
        // Get the linear solver
        libMesh::LinearSolver<libMesh::Number> *linear_solver = context.system->get_linear_solver();

        // Set ourselves to reuse the preconditioner
        linear_solver->reuse_preconditioner(true);

        // Solve the adjoint problem
        this->steady_adjoint_solve(context);

        // Go back to not reusing the preconditioner
        linear_solver->reuse_preconditioner(false);
      }

    if( context.output_adjoint )
      context.vis->output_adjoint( context.equation_system, context.system );

    if( context.output_vis )
      {
        context.postprocessing->update_quantities( *(context.equation_system) );
        context.vis->output( context.equation_system );
      }

    if( context.output_residual ) context.vis->output_residual( context.equation_system, context.system );

    return;
  }

  void SteadySolver::adjoint_qoi_parameter_sensitivity
  (SolverContext& context,
   const libMesh::QoISet&          qoi_indices,
   const libMesh::ParameterVector& parameters_in,
   libMesh::SensitivityData&       sensitivities) const
  {
    // Get the linear solver
    libMesh::LinearSolver<libMesh::Number> *linear_solver = context.system->get_linear_solver();

    // Set ourselves to reuse the preconditioner
    linear_solver->reuse_preconditioner(true);

    context.system->adjoint_qoi_parameter_sensitivity
      (qoi_indices, parameters_in, sensitivities);

    // Go back to not reusing the preconditioner
    linear_solver->reuse_preconditioner(false);
  }

  void SteadySolver::forward_qoi_parameter_sensitivity
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

} // namespace GRINS
