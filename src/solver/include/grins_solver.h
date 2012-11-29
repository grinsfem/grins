//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
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
#include "multiphysics_sys.h"
#include "visualization.h"
#include "bc_factory.h"

// libMesh
#include "getpot.h"
#include "libmesh.h"
#include "libmesh_logging.h"
#include "mesh.h"
#include "diff_solver.h"
#include "equation_systems.h"
#include "euler_solver.h"
#include "steady_solver.h"
#include "boundary_conditions.h"

#ifdef HAVE_GRVY
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
    
    virtual void solve( GRINS::MultiphysicsSystem* system,
			std::tr1::shared_ptr<libMesh::EquationSystems> equation_system =
			std::tr1::shared_ptr<libMesh::EquationSystems>(),
			std::tr1::shared_ptr<GRINS::Visualization> vis = 
			std::tr1::shared_ptr<GRINS::Visualization>(),
			bool output_vis = false,
			bool output_residual = false )=0;

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
