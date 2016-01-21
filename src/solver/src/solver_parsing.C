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

// These functions
#include "grins/solver_names.h"
#include "grins/solver_parsing.h"
#include "grins/strategies_parsing.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  std::string SolverParsing::solver_type(const GetPot& input)
  {
    bool mesh_adaptive = StrategiesParsing::is_mesh_adaptive(input);

    bool transient = SolverParsing::is_transient(input);

    std::string solver_type = input("SolverOptions/solver_type", "DIE!");

    if( solver_type == SolverNames::displacement_continuation() )
      {
        // Don't need to do anything
      }
    else if(transient && !mesh_adaptive)
      {
        solver_type = SolverNames::unsteady_solver();
      }
    else if( !transient && !mesh_adaptive )
      {
        solver_type = SolverNames::steady_solver();
      }
    else if( !transient && mesh_adaptive )
      {
        solver_type = SolverNames::steady_mesh_adaptive_solver();
      }
    else if( transient && mesh_adaptive )
      {
        solver_type = SolverNames::unsteady_mesh_adaptive_solver();
      }
    else
      {
        libmesh_error_msg("Unsupported combination of solver options!");
      }

    return solver_type;
  }

  bool SolverParsing::is_transient( const GetPot& input )
  {
    return input("unsteady-solver/transient", false );
  }

} // end namespace GRINS
