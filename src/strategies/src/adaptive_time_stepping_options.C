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

// This class
#include "grins/adaptive_time_stepping_options.h"

// GRINS
#include "grins/common.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"

namespace GRINS
{
  AdaptiveTimeSteppingOptions::AdaptiveTimeSteppingOptions( const GetPot& input )
    : _is_time_adaptive(false),
      _target_tolerance(0.0),
      _upper_tolerance(0.0),
      _max_growth(0.0)
  {
    this->check_dup_input_style(input);

    if( this->is_old_style(input) )
      this->parse_old_style(input);
    else
      this->parse_new_style(input);

    if( _target_tolerance != 0.0 )
      _is_time_adaptive = true;
  }

  void AdaptiveTimeSteppingOptions::check_dup_input_style( const GetPot& input ) const
  {
    if( (input.have_variable("unsteady-solver/target_tolerance") &&
         input.have_section("Strategies/AdaptiveTimeStepping/target_tolerance")) ||
        (input.have_variable("unsteady-solver/upper_tolerance") &&
         input.have_section("Strategies/AdaptiveTimeStepping/upper_tolerance")) ||
        (input.have_variable("unsteady-solver/max_growth") &&
         input.have_section("Strategies/AdaptiveTimeStepping/max_growth")) )
      libmesh_error_msg("ERROR: Cannot use both old and new style of options for AdaptiveTimeSteppingOptions!");
  }

  bool AdaptiveTimeSteppingOptions::is_old_style( const GetPot& input ) const
  {
    return input.have_variable("unsteady-solver/target_tolerance");
  }

  void AdaptiveTimeSteppingOptions::parse_old_style(const GetPot& input)
  {
    {
      std::string warning = "WARNING: Using [MeshAdaptivity/<options>] is a DEPRECATED\n";
      warning += "         style of input for ErrorEstimator options. Please\n";
      warning += "         update to use the [Strategies/ErrorEstimation/<options> style.\n";
      grins_warning(warning);
    }

    std::string section = "unsteady-solver";
    this->parse_options(input,section);
  }

  void AdaptiveTimeSteppingOptions::parse_new_style(const GetPot& input)
  {
    std::string section = "Strategies/AdaptiveTimeStepping";
    this->parse_options(input,section);
  }

  void AdaptiveTimeSteppingOptions::parse_options(const GetPot& input, const std::string& section)
  {
    // If the user set the target tolerance, for them to set the other values too
    if( input.have_variable(section+"/target_tolerance") )
      {
        if( !input.have_variable(section+"/upper_tolerance") )
          libmesh_error_msg("ERROR: Must specify "+section+"/upper_tolerance for adaptive time stepping!");

        if( !input.have_variable(section+"/max_growth") )
          libmesh_error_msg("ERROR: Must specify "+section+"/max_growth for adaptive time stepping!");
      }

    _target_tolerance = input(section+"/target_tolerance", 0.0 );
    _upper_tolerance = input(section+"/upper_tolerance", 0.0 );
    _max_growth = input(section+"/max_growth", 0.0 );

    // parse component_norm
    const unsigned int n_component_norm =
      input.vector_variable_size(section+"/component_norm");

    for (unsigned int i=0; i != n_component_norm; ++i)
      {
        const std::string current_norm = input(section+"/component_norm", std::string("L2"), i);
        _component_norm.set_type(i, libMesh::Utility::string_to_enum<libMesh::FEMNormType>(current_norm) );
      }
  }

} // end namespace GRINS
