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
#include "grins/time_stepping_parsing.h"

// GRINS
#include "grins/common.h"
#include "grins/solver_parsing.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  unsigned int TimeSteppingParsing::parse_n_timesteps( const GetPot& input )
  {
    SolverParsing::dup_solver_option_check(input,
                                           "unsteady-solver/n_timesteps",
                                           "SolverOptions/TimeStepping/n_timesteps");

    unsigned int n_timesteps = 0;

    if( input.have_variable("unsteady-solver/n_timesteps") )
      {
        n_timesteps = input("unsteady-solver/n_timesteps",0);

        std::string warning = "WARNING: unsteady-solver/n_timesteps is DEPRECATED!\n";
        warning += "        Please use SolverOptions/TimeStepping/n_timesteps to specify # of timesteps.\n";
        grins_warning(warning);
      }
    else if( input.have_variable("SolverOptions/TimeStepping/n_timesteps") )
        n_timesteps = input("SolverOptions/TimeStepping/n_timesteps",0);
    else
      libmesh_error_msg("ERROR: Could not find valid entry for n_timesteps!");

    return n_timesteps;
  }

  unsigned int TimeSteppingParsing::parse_backtrack_deltat( const GetPot& input )
  {
    SolverParsing::dup_solver_option_check(input,
                                           "unsteady-solver/backtrack_deltat",
                                           "SolverOptions/TimeStepping/backtrack_deltat");

    unsigned int backtrack_deltat = 0;

    if( input.have_variable("unsteady-solver/backtrack_deltat") )
      {
        backtrack_deltat = input("unsteady-solver/backtrack_deltat",0);

        std::string warning = "WARNING: unsteady-solver/backtrack_deltat is DEPRECATED!\n";
        warning += "        Please use SolverOptions/TimeStepping/backtrack_deltat to set backtrack_deltat.\n";
        grins_warning(warning);
      }
    else
      backtrack_deltat = input("SolverOptions/TimeStepping/backtrack_deltat",0);

    return backtrack_deltat;
  }

  double TimeSteppingParsing::parse_theta( const GetPot& input )
  {
    return input("unsteady-solver/theta", 0.5 );
  }

  double TimeSteppingParsing::parse_deltat( const GetPot& input )
  {
    return input("unsteady-solver/deltat", 0.0 );
  }

} // end namespace GRINS
