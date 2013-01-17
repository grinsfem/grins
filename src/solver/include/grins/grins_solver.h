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

#ifndef GRINS_SOLVER_H
#define GRINS_SOLVER_H

#include "boost/tr1/memory.hpp"

// GRINS
#include "grins_config.h"
#include "grins/multiphysics_sys.h"
#include "grins/visualization.h"
#include "grins/bc_factory.h"
#include "grins/solver_context.h"

// libMesh
#include "getpot.h"
#include "libmesh.h"
#include "libmesh_logging.h"
#include "mesh.h"
#include "diff_solver.h"
#include "equation_systems.h"
#include "euler_solver.h"
#include "steady_solver.h"
#include "grins/boundary_conditions.h"

#ifdef GRINS_HAVE_GRVY
#include "grvy.h" // GRVY timers
#endif

namespace GRINS
{
  class Solver
  {
  public:
    Solver( const GetPot& input );
    virtual ~Solver();

    virtual void initialize( const GetPot& input, 
			     std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
			     GRINS::MultiphysicsSystem* system );
    
    virtual void solve( SolverContext& context )=0;

  protected:

    // Linear/Nonlinear solver options
    unsigned int _max_nonlinear_iterations;
    double _relative_step_tolerance;
    double _absolute_step_tolerance;
    double _relative_residual_tolerance;
    double _absolute_residual_tolerance;
    unsigned int _max_linear_iterations;
    double _initial_linear_tolerance;

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
