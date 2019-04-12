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

#ifndef GRINS_TIME_STEPPING_PARSING_H
#define GRINS_TIME_STEPPING_PARSING_H

// C++
#include <string>

// Forward declarations
class GetPot;

namespace GRINS
{
  class TimeSteppingParsing
  {
  public:

    TimeSteppingParsing(){};

    ~TimeSteppingParsing(){};

    static unsigned int parse_n_timesteps( const GetPot& input );

    //! Parse option to retry failed time steps with smaller \f$ \Delta t \f$
    /*! backtrack_deltat is the number of time the TimeSolver will try
      to resolve the timestep with a smaller \f$ \Delta t \f$. Default is 0. */
    static unsigned int parse_backtrack_deltat( const GetPot& input );

    //! Parse value of \f$ \theta \f$ for theta method time stepping.
    /*! \f$ \theta \in [0,1] \f$ Option only used for theta method-based
      time solvers. Defaults to 0.5. */
    static double parse_theta( const GetPot& input );

    static double parse_deltat( const GetPot& input );

    static std::string parse_time_stepper_name( const GetPot& input );

    static bool parse_zero_initial_velocity( const GetPot & input );

    static bool parse_recompute_accel( const GetPot & input );
  };

} // end namespace GRINS

#endif // GRINS_TIME_STEPPING_PARSING_H
