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

#ifndef GRINS_MESH_ADAPTIVE_SOLVER_BASE_H
#define GRINS_MESH_ADAPTIVE_SOLVER_BASE_H

// C++
#include <string>

//libMesh
#include "libmesh/libmesh.h"
#include "libmesh/mesh_refinement.h"

// libMesh forward declarations
class GetPot;
namespace libMesh
{
  class MeshBase;
  class ErrorVector;
}

namespace GRINS
{
  class SolverContext;

  class MeshAdaptiveSolverBase
  {
  public:

    MeshAdaptiveSolverBase( const GetPot& input );

    virtual ~MeshAdaptiveSolverBase();

  protected:

    enum RefinementFlaggingType{ INVALID = 0,
                                 ERROR_TOLERANCE,
                                 N_ELEM_TARGET,
                                 ERROR_FRACTION,
                                 ELEM_FRACTION,
                                 MEAN_STD_DEV };

    unsigned int _max_refinement_steps;
    bool _coarsen_by_parents;
    libMesh::Real _absolute_global_tolerance;
    unsigned int _nelem_target;
    libMesh::Real _refine_fraction;
    libMesh::Real _coarsen_fraction;
    libMesh::Real _coarsen_threshold;
    bool _compute_qoi_error_estimate;
    bool _plot_cell_errors;
    std::string _error_plot_prefix;

    unsigned int _node_level_mismatch_limit;
    unsigned int _edge_level_mismatch_limit;
    unsigned int _face_level_mismatch_limit;
    bool _enforce_mismatch_limit_prior_to_refinement;

    unsigned int _max_h_level;

    RefinementFlaggingType _refinement_type;

    libMesh::UniquePtr<libMesh::MeshRefinement> _mesh_refinement;

    void build_mesh_refinement( libMesh::MeshBase& mesh );

    void set_refinement_type( const GetPot& input,
                              RefinementFlaggingType& refinement_type );

    bool check_for_convergence( SolverContext& context,
                                const libMesh::ErrorVector& error ) const;

    void flag_elements_for_refinement( const libMesh::ErrorVector& error );

    void estimate_error_for_amr( SolverContext& context, libMesh::ErrorVector& error );

    void perform_amr( SolverContext& context, const libMesh::ErrorVector& error );

  private:

    MeshAdaptiveSolverBase();

  };

} // end namespace GRINS

#endif // GRINS_MESH_ADAPTIVE_SOLVER_BASE_H
