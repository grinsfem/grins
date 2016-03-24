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

#ifndef GRINS_MESH_ADAPTIVITY_OPTIONS_H
#define GRINS_MESH_ADAPTIVITY_OPTIONS_H

// C++
#include <string>

// libMesh
#include "libmesh/libmesh.h"

// libmMesh forward declarations
class GetPot;

namespace GRINS
{
  //! Container for mesh adaptivity options
  class MeshAdaptivityOptions
  {
  public:
    MeshAdaptivityOptions( const GetPot& input );
    ~MeshAdaptivityOptions(){};

    bool is_mesh_adaptive() const
    { return _is_mesh_adaptive; }

    const std::string& refinement_strategy() const
    { return _refinement_strategy; }

    unsigned int max_refinement_steps() const
    { return _max_refinement_steps; }

    bool coarsen_by_parents() const
    { return _coarsen_by_parents; }

    libMesh::Real absolute_global_tolerance() const
    { return _absolute_global_tolerance; }

    unsigned int nelem_target() const
    { return _nelem_target; }

    libMesh::Real refine_fraction() const
    { return _refine_fraction; }

    libMesh::Real coarsen_fraction () const
    { return _coarsen_fraction; }

    libMesh::Real coarsen_threshold() const
    { return _coarsen_threshold; }

    bool plot_cell_errors() const
    { return _plot_cell_errors; }

    const std::string& error_plot_prefix() const
    { return _error_plot_prefix; }

    unsigned int node_level_mismatch_limit() const
    { return _node_level_mismatch_limit; }

    unsigned int edge_level_mismatch_limit() const
    { return _edge_level_mismatch_limit; }

    unsigned int face_level_mismatch_limit() const
    { return _face_level_mismatch_limit; }

    bool enforce_mismatch_limit_prior_to_refinement() const
    { return _enforce_mismatch_limit_prior_to_refinement; }

    unsigned int max_h_level() const
    { return _max_h_level; }

  private:

    void check_dup_input_style( const GetPot& input ) const;

    bool is_old_style( const GetPot& input ) const;

    void parse_old_style(const GetPot& input);

    void parse_new_style(const GetPot& input);

    void parse_options(const GetPot& input, const std::string& section);

    bool   _is_mesh_adaptive;
    std::string _refinement_strategy;
    unsigned int _max_refinement_steps;
    bool _coarsen_by_parents;
    libMesh::Real _absolute_global_tolerance;
    unsigned int _nelem_target;
    libMesh::Real _refine_fraction;
    libMesh::Real _coarsen_fraction;
    libMesh::Real _coarsen_threshold;
    bool _plot_cell_errors;
    std::string _error_plot_prefix;

    unsigned int _node_level_mismatch_limit;
    unsigned int _edge_level_mismatch_limit;
    unsigned int _face_level_mismatch_limit;
    bool _enforce_mismatch_limit_prior_to_refinement;

    unsigned int _max_h_level;

  };

} // end namespace GRINS

#endif // GRINS_MESH_ADAPTIVITY_OPTIONS_H
