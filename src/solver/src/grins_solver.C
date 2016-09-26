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


// C++
#include <iostream>
#include <iomanip>

// This class
#include "grins/grins_solver.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/solver_context.h"
#include "grins/composite_qoi.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/diff_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/dof_map.h"

namespace GRINS
{

  Solver::Solver( const GetPot& input )
    : _max_nonlinear_iterations( input("linear-nonlinear-solver/max_nonlinear_iterations", 10 ) ),
      _relative_step_tolerance( input("linear-nonlinear-solver/relative_step_tolerance", 1.e-6 ) ),
      _absolute_step_tolerance( input("linear-nonlinear-solver/absolute_step_tolerance", 0.0 ) ),
      _relative_residual_tolerance( input("linear-nonlinear-solver/relative_residual_tolerance", 1.e-15 ) ),
      _absolute_residual_tolerance( input("linear-nonlinear-solver/absolute_residual_tolerance", 0.0 ) ),
      _initial_linear_tolerance( input("linear-nonlinear-solver/initial_linear_tolerance", 1.e-3 ) ),
      _minimum_linear_tolerance( input("linear-nonlinear-solver/minimum_linear_tolerance", 1.e-3 ) ),
      _max_linear_iterations( input("linear-nonlinear-solver/max_linear_iterations", 500 ) ),
      _continue_after_backtrack_failure( input("linear-nonlinear-solver/continue_after_backtrack_failure", false ) ),
      _continue_after_max_iterations( input("linear-nonlinear-solver/continue_after_max_iterations", false ) ),
      _require_residual_reduction( input("linear-nonlinear-solver/require_residual_reduction", true ) ),
      _solver_quiet( input("screen-options/solver_quiet", false ) ),
      _solver_verbose( input("screen-options/solver_verbose", false ) )
  {
    return;
  }


  Solver::~Solver()
  {
    return;
  }

  void Solver::initialize( const GetPot& /*input*/,
			   SharedPtr<libMesh::EquationSystems> equation_system,
			   MultiphysicsSystem* system )
  {

    // Defined in subclasses depending on the solver used.
    this->init_time_solver(system);

    // Initialize the system
    equation_system->init();

    // Get diff solver to set options
    libMesh::DiffSolver &solver = *(system->time_solver->diff_solver().get());

    // Set linear/nonlinear solver options
    this->set_solver_options( solver );

    return;
  }

  void Solver::set_solver_options( libMesh::DiffSolver& solver  )
  {
    solver.quiet                       = this->_solver_quiet;
    solver.verbose                     = this->_solver_verbose;
    solver.max_nonlinear_iterations    = this->_max_nonlinear_iterations;
    solver.relative_step_tolerance     = this->_relative_step_tolerance;
    solver.absolute_step_tolerance     = this->_absolute_step_tolerance;
    solver.relative_residual_tolerance = this->_relative_residual_tolerance;
    solver.absolute_residual_tolerance = this->_absolute_residual_tolerance;
    solver.max_linear_iterations       = this->_max_linear_iterations;
    solver.continue_after_backtrack_failure = this->_continue_after_backtrack_failure;
    solver.initial_linear_tolerance    = this->_initial_linear_tolerance;
    solver.minimum_linear_tolerance    = this->_minimum_linear_tolerance;
    solver.continue_after_max_iterations    = this->_continue_after_max_iterations;
    if(dynamic_cast<libMesh::NewtonSolver*>(&solver))
      {
        dynamic_cast<libMesh::NewtonSolver&>(solver).require_residual_reduction = this->_require_residual_reduction;
      }
    else
      {
        // If the user tried to set require_residual_reduction flag to false
        // despite not having a NewtonSolver spit out a warning
        if(this->_require_residual_reduction == false)
          {
            libmesh_warning("GRINS can't change require_residual_reduction when not using NewtonSolver!");
          }
      }

    return;
  }

  void Solver::steady_adjoint_solve( SolverContext& context )
  {
    libMesh::out << "==========================================================" << std::endl
                 << "Solving adjoint problem." << std::endl
                 << "==========================================================" << std::endl;

    context.system->adjoint_solve();
    context.system->set_adjoint_already_solved(true);
  }

  void Solver::print_scalar_vars( SolverContext& context )
  {
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
  }

  void Solver::print_qoi( SolverContext & context )
  {
    context.system->assemble_qoi();
    const CompositeQoI* my_qoi = libMesh::cast_ptr<const CompositeQoI*>(context.system->get_qoi());
    context.qoi_output->output_qois(*my_qoi, context.system->comm());
  }

} // namespace GRINS
