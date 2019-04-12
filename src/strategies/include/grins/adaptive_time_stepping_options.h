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

#ifndef GRINS_ADAPTIVE_TIME_STEPPING_OPTIONS_H
#define GRINS_ADAPTIVE_TIME_STEPPING_OPTIONS_H

// C++
#include <string>

// libMesh
#include "libmesh/system_norm.h"

// libmMesh forward declarations
class GetPot;

namespace GRINS
{
  //! Container for adaptive time-stepping options
  class AdaptiveTimeSteppingOptions
  {
  public:
    AdaptiveTimeSteppingOptions( const GetPot& input );
    ~AdaptiveTimeSteppingOptions(){};

    bool is_time_adaptive() const
    { return _is_time_adaptive; }

    double target_tolerance() const
    { return _target_tolerance; }

    double upper_tolerance() const
    { return _upper_tolerance; }

    double max_growth() const
    { return _max_growth; }

    const libMesh::SystemNorm& component_norm()
    { return _component_norm; }

  private:

    void check_dup_input_style( const GetPot& input ) const;

    bool is_old_style( const GetPot& input ) const;

    void parse_old_style(const GetPot& input);

    void parse_new_style(const GetPot& input);

    void parse_options(const GetPot& input, const std::string& section);

    bool   _is_time_adaptive;

    // target tolerance parameter for adaptive time stepping
    /*! 0.0 means there is no adaptive time stepping enabled. To enable
      adaptive time stepping with the libMesh::TwostepTimeSolver, this
      parameter should be positive. */
    double _target_tolerance;
    double _upper_tolerance;
    double _max_growth;

    // The norm to use for each solution variable for adaptive time stepping
    libMesh::SystemNorm _component_norm;
  };

} // end namespace GRINS

#endif // GRINS_ADAPTIVE_TIME_STEPPING_OPTIONS_H
