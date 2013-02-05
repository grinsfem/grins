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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// C++
#include <iostream>

// This class
#include "grins/grins_solver.h"

// GRINS
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/diff_solver.h"

namespace GRINS
{

  Solver::Solver( const GetPot& input )
    : _max_nonlinear_iterations( input("linear-nonlinear-solver/max_nonlinear_iterations", 10 ) ),
      _relative_step_tolerance( input("linear-nonlinear-solver/relative_step_tolerance", 1.e-6 ) ),
      _absolute_step_tolerance( input("linear-nonlinear-solver/absolute_step_tolerance", 0.0 ) ),
      _relative_residual_tolerance( input("linear-nonlinear-solver/relative_residual_tolerance", 1.e-15 ) ),
      _absolute_residual_tolerance( input("linear-nonlinear-solver/absolute_residual_tolerance", 0.0 ) ),
      _max_linear_iterations( input("linear-nonlinear-solver/max_linear_iterations", 500 ) ),
      _initial_linear_tolerance( input("linear-nonlinear-solver/initial_linear_tolerance", 1.e-3 ) ),
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
			   std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
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
    solver.initial_linear_tolerance    = this->_initial_linear_tolerance;

    return;
  }

} // namespace GRINS
