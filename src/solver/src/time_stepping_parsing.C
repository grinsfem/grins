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

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  unsigned int TimeSteppingParsing::parse_n_timesteps( const GetPot& input )
  {
    return input("unsteady-solver/n_timesteps", 1 );
  }

  unsigned int TimeSteppingParsing::parse_backtrack_deltat( const GetPot& input )
  {
    return input("unsteady-solver/backtrack_deltat", 0 );
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
