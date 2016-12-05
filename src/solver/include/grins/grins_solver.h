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


#ifndef GRINS_SOLVER_H
#define GRINS_SOLVER_H

// GRINS
#include "grins/shared_ptr.h"

// libMesh
#include "libmesh/equation_systems.h"

#ifdef GRINS_HAVE_GRVY
#include "libmesh/ignore_warnings.h" // avoid auto_ptr deprecated warnings
#include "grvy.h" // GRVY timers
#include "libmesh/restore_warnings.h"
#endif

// libMesh forward declarations
class GetPot;

namespace libMesh
{
  class DiffSolver;
  class ParameterVector;
  class SensitivityData;
}

namespace GRINS
{
  // Forward declarations
  class MultiphysicsSystem;
  class SolverContext;

  class Solver
  {
  public:
    Solver( const GetPot& input );
    virtual ~Solver();

    virtual void initialize( const GetPot& input,
			     SharedPtr<libMesh::EquationSystems> equation_system,
			     GRINS::MultiphysicsSystem* system );

    virtual void solve( SolverContext& context )=0;

    virtual void adjoint_qoi_parameter_sensitivity
      (SolverContext&                  /*context*/,
       const libMesh::QoISet&          /*qoi_indices*/,
       const libMesh::ParameterVector& /*parameters_in*/,
       libMesh::SensitivityData&       /*sensitivities*/)
      const
    { libmesh_not_implemented(); }

    virtual void forward_qoi_parameter_sensitivity
      (SolverContext&                  /*context*/,
       const libMesh::QoISet&          /*qoi_indices*/,
       const libMesh::ParameterVector& /*parameters_in*/,
       libMesh::SensitivityData&       /*sensitivities*/)
      const
    { libmesh_not_implemented(); }

    //! Do steady version of adjoint solve
    /*! We put this here since we may want to reuse this
        in multiple different steady solves. */
    void steady_adjoint_solve( SolverContext& context );

    void print_scalar_vars( SolverContext& context );

    void print_qoi( SolverContext& context );

  protected:

    // Linear/Nonlinear solver options
    unsigned int _max_nonlinear_iterations;
    double _relative_step_tolerance;
    double _absolute_step_tolerance;

    // _relative_residual_tolerance applies to both one of the
    // stopping criteria for (nonlinear) forward solves and *the*
    // stopping criterion for (linear) adjoint solves.
    double _relative_residual_tolerance;

    double _absolute_residual_tolerance;
    double _initial_linear_tolerance;
    double _minimum_linear_tolerance;
    unsigned int _max_linear_iterations;
    bool _continue_after_backtrack_failure;
    bool _continue_after_max_iterations;
    bool _require_residual_reduction;

    // Screen display options
    bool _solver_quiet;
    bool _solver_verbose;

    void set_solver_options( libMesh::DiffSolver& solver );

    virtual void init_time_solver(GRINS::MultiphysicsSystem* system)=0;

  };

} //End namespace block

#endif //GRINS_SOLVER_H
