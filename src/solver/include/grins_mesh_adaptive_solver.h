//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_MESH_ADAPTIVE_SOLVER_H
#define GRINS_MESH_ADAPTIVE_SOLVER_H

//GRINS
#include "grins_solver.h"
#include "qoi_base.h"
#include "visualization.h"

//libMesh
#include "auto_ptr.h"
#include "equation_systems.h"
#include "error_vector.h"
#include "mesh_refinement.h"

namespace GRINS
{
  class MeshAdaptiveSolver : public Solver
  {
  public:
    MeshAdaptiveSolver( const GetPot& input );

    virtual ~MeshAdaptiveSolver();

    virtual void solve( GRINS::MultiphysicsSystem* system,
			std::tr1::shared_ptr<libMesh::EquationSystems> equation_system =
			std::tr1::shared_ptr<libMesh::EquationSystems>(),
      std::tr1::shared_ptr<GRINS::QoIBase> qoi_base =
      std::tr1::shared_ptr<GRINS::QoIBase>(),
			std::tr1::shared_ptr<GRINS::Visualization> vis = 
			std::tr1::shared_ptr<GRINS::Visualization>(),
			bool output_vis = false,
			bool output_residual = false,
      std::tr1::shared_ptr<libMesh::ErrorEstimator> =
      std::tr1::shared_ptr<libMesh::ErrorEstimator>() );

  protected:
    unsigned int _max_r_steps;
    bool _coarsen_by_parents;
    Real _absolute_global_tolerance;
    unsigned int _nelem_target;
    Real _refine_fraction;
    Real _coarsen_fraction;
    Real _coarsen_threshold;
    bool _output_adjoint_sol;
    bool _plot_cell_errors;
    std::string _error_plot_prefix;

    libMesh::AutoPtr<MeshRefinement> _mesh_refinement;

    virtual void init_time_solver( GRINS::MultiphysicsSystem* system );

    void read_input_options( const GetPot& input );

    void build_mesh_refinement( MeshBase &mesh );
  };
} // namespace GRINS
#endif // GRINS_STEADY_SOLVER_H
