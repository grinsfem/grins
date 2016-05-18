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


#ifndef GRINS_SIMULATION_BUILDER_H
#define GRINS_SIMULATION_BUILDER_H

// GRINS
#include "grins/mesh_builder.h"
#include "grins/solver_factory.h"
#include "grins/visualization_factory.h"
#include "grins/qoi_factory.h"
#include "grins/postprocessing_factory.h"
#include "grins/shared_ptr.h"

namespace GRINS
{
  class SimulationBuilder
  {
  public:

    SimulationBuilder();
    virtual ~SimulationBuilder(){};

    SharedPtr<libMesh::UnstructuredMesh> build_mesh
      ( const GetPot& input,
        const libMesh::Parallel::Communicator &comm
        LIBMESH_CAN_DEFAULT_TO_COMMWORLD );

    SharedPtr<GRINS::Solver> build_solver( const GetPot& input );

    SharedPtr<GRINS::Visualization> build_vis
      ( const GetPot& input,
        const libMesh::Parallel::Communicator &comm
        LIBMESH_CAN_DEFAULT_TO_COMMWORLD );

    SharedPtr<CompositeQoI> build_qoi( const GetPot& input );

    SharedPtr<PostProcessedQuantities<libMesh::Real> > build_postprocessing( const GetPot& input );

    void attach_solver_factory( SharedPtr<SolverFactory> solver_factory );

    void attach_mesh_builder( SharedPtr<MeshBuilder> mesh_builder );

    void attach_vis_factory( SharedPtr<VisualizationFactory> vis_factory );

    void attach_qoi_factory( SharedPtr<QoIFactory> qoi_factory );

    void attach_postprocessing_factory( SharedPtr<PostprocessingFactory> postprocessing_factory );

    const MeshBuilder& mesh_builder() const;

  protected:

    SharedPtr<MeshBuilder> _mesh_builder;
    SharedPtr<SolverFactory> _solver_factory;
    SharedPtr<VisualizationFactory> _vis_factory;
    SharedPtr<QoIFactory> _qoi_factory;
    SharedPtr<PostprocessingFactory> _postprocessing_factory;

  }; //class SimulationBuilder
} // namespace GRINS

#endif //GRINS_SIMULATION_BUILDER_H
