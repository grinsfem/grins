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

// These functions
#include "grins/solver_names.h"
#include "grins/solver_parsing.h"
#include "grins/mesh_adaptivity_options.h"

// GRINS
#include "grins/common.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  std::string SolverParsing::solver_type(const GetPot& input)
  {
    std::string solver_type;

    // If the user set the solver_type, use that
    if( input.have_variable("SolverOptions/solver_type") )
      solver_type = input("SolverOptions/solver_type", "DIE!");

    // Otherwise, they're using grins_steady/unsteady (possibly mesh_adaptive) solver
    else
      {
        //! \todo We should just pass this object in from the calling function
        MeshAdaptivityOptions mesh_adaptivity_options(input);
        bool mesh_adaptive = mesh_adaptivity_options.is_mesh_adaptive();

        bool transient = SolverParsing::is_transient(input);

        if(transient && !mesh_adaptive)
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
      }

    return solver_type;
  }

  bool SolverParsing::is_transient( const GetPot& input )
  {
    // Can't specify both old and new version
    SolverParsing::dup_solver_option_check(input,
                                           "unsteady-solver/transient",
                                           "SolverOptions/TimeStepping/solver_type");

    bool transient = false;

    if( input.have_variable("unsteady-solver/transient") )
      {
        transient = input("unsteady-solver/transient",false);

        std::string warning = "WARNING: unsteady-solver/transient is DEPRECATED!\n";
        warning += "        Please use SolverOptions/TimeStepping/solver_type to specify time stepping solver.\n";
        grins_warning(warning);
      }

    // In the new version, we set the solver type so we just need to
    // check if the variable is present. Solver class will figure
    // out the type.
    if( input.have_variable("SolverOptions/TimeStepping/solver_type") )
      transient = true;

    return transient;
  }

  void SolverParsing::dup_solver_option_check( const GetPot& input,
                                               const std::string& option1,
                                               const std::string& option2 )
  {
    // Can't specify both old and new version
    if( input.have_variable(option1) && input.have_variable(option2) )
      {
        libmesh_error_msg("ERROR: Cannot specify both "+option1+" and "+option2);
      }
  }

} // end namespace GRINS
