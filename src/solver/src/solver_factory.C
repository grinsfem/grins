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
#include "grins/solver_factory.h"

// GRINS
#include "grins/solver_names.h"
#include "grins/solver_parsing.h"
#include "grins/grins_steady_solver.h"
#include "grins/grins_unsteady_solver.h"
#include "grins/steady_mesh_adaptive_solver.h"
#include "grins/unsteady_mesh_adaptive_solver.h"
#include "grins/displacement_continuation_solver.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  SharedPtr<Solver> SolverFactory::build(const GetPot& input)
  {
    std::string solver_type = SolverParsing::solver_type(input);

    SharedPtr<Solver> solver;  // Effectively NULL

    if( solver_type == SolverNames::displacement_continuation() )
      {
        solver.reset( new DisplacementContinuationSolver(input) );
      }
    else if(solver_type == SolverNames::unsteady_solver() )
      {
        solver.reset( new UnsteadySolver(input) );
      }
    else if( solver_type == SolverNames::steady_solver() )
      {
        solver.reset( new SteadySolver(input) );
      }
    else if( solver_type == SolverNames::steady_mesh_adaptive_solver() )
      {
        solver.reset( new SteadyMeshAdaptiveSolver(input) );
      }
    else if( solver_type == SolverNames::unsteady_mesh_adaptive_solver() )
      {
        solver.reset( new UnsteadyMeshAdaptiveSolver(input) );
      }
    else
      {
        libmesh_error_msg("Invalid solver_type: "+solver_type);
      }

    return solver;
  }

} // namespace GRINS
