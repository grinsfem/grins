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
#include "grins/mesh_adaptivity_options.h"

// GRINS
#include "grins/common.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"

namespace GRINS
{
  MeshAdaptivityOptions::MeshAdaptivityOptions( const GetPot& input )
    : _is_mesh_adaptive(false),
      _refinement_strategy("elem_fraction"),
      _max_refinement_steps(0),
      _coarsen_by_parents(true),
      _absolute_global_tolerance(0),
      _nelem_target(0),
      _refine_fraction(0.2),
      _coarsen_fraction(0.2),
      _coarsen_threshold(0),
      _plot_cell_errors(false),
      _error_plot_prefix("cell_error"),
      _node_level_mismatch_limit(0),
      _edge_level_mismatch_limit(0),
      _face_level_mismatch_limit(1),
      _enforce_mismatch_limit_prior_to_refinement(false),
      _max_h_level(libMesh::invalid_uint)
  {
    this->check_dup_input_style(input);

    if( this->is_old_style(input) )
      this->parse_old_style(input);
    else
      this->parse_new_style(input);
  }

  void MeshAdaptivityOptions::check_dup_input_style( const GetPot& input ) const
  {
    if( input.have_variable("MeshAdaptivity") &&
        input.have_section("Strategies/MeshAdaptivity") )
      libmesh_error_msg("ERROR: Cannot use both old and new style of options for MeshAdaptivityOptions!");
  }

  bool MeshAdaptivityOptions::is_old_style( const GetPot& input ) const
  {
    return input.have_section("MeshAdaptivity");
  }

  void MeshAdaptivityOptions::parse_old_style(const GetPot& input)
  {
    {
      std::string warning = "WARNING: Using [MeshAdaptivity/<options>] is a DEPRECATED\n";
      warning += "         style of input for mesh adaptivity options. Please\n";
      warning += "         update to use the [Strategies/MeshAdaptivity/<options> style.\n";
      grins_warning(warning);
    }

    std::string section = "MeshAdaptivity";
    this->parse_options(input,section);
  }

  void MeshAdaptivityOptions::parse_new_style(const GetPot& input)
  {
    std::string section = "Strategies/MeshAdaptivity";
    this->parse_options(input,section);
  }

  void MeshAdaptivityOptions::parse_options(const GetPot& input, const std::string& section)
  {
    _is_mesh_adaptive = input(section+"/mesh_adaptive", false);
    _refinement_strategy = input(section+"/refinement_strategy", "elem_fraction" );
    _max_refinement_steps = input(section+"/max_refinement_steps", 0);
    _coarsen_by_parents = true; //! \todo PB: why aren't parsing this?
    _absolute_global_tolerance = input(section+"/absolute_global_tolerance", 0);
    _nelem_target = input(section+"/nelem_target", 0);
    _refine_fraction = input(section+"/refine_percentage", 0.8);
    _coarsen_fraction = input(section+"/coarsen_percentage", 0.1);
    _coarsen_threshold = input(section+"/coarsen_threshold", 0);
    _plot_cell_errors = input(section+"/plot_cell_errors", false);
    _error_plot_prefix = input(section+"/error_plot_prefix", "cell_error");
    _node_level_mismatch_limit = input(section+"/node_level_mismatch_limit", 0);
    _edge_level_mismatch_limit = input(section+"/edge_level_mismatch_limit", 0 );
    _face_level_mismatch_limit = input(section+"/face_level_mismatch_limit", 1 );
    _enforce_mismatch_limit_prior_to_refinement = input(section+"/enforce_mismatch_limit_prior_to_refinement", true );
    _max_h_level = input(section+"/max_h_level",libMesh::invalid_uint);
  }

} // end namespace GRINS
