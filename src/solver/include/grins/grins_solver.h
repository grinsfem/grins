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


#ifndef GRINS_SOLVER_H
#define GRINS_SOLVER_H

// C++
#include "boost/tr1/memory.hpp"

// GRINS
#include "grins/nbc_container.h"

// libMesh
#include "libmesh/equation_systems.h"

#ifdef GRINS_HAVE_GRVY
#include "grvy.h" // GRVY timers
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
			     std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
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

    // Screen display options
    bool _solver_quiet;
    bool _solver_verbose;    
    
    /* Keep copies of the boundary conditions around
       in case they need to be updated during a solve;
       for example parameter continuation. */
    std::map< std::string, GRINS::NBCContainer > _neumann_bc_funcs;

    void set_solver_options( libMesh::DiffSolver& solver );

    virtual void init_time_solver(GRINS::MultiphysicsSystem* system)=0;

  };

} //End namespace block

#endif //GRINS_SOLVER_H
