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


#ifndef GRINS_SIMULATION_BUILDER_H
#define GRINS_SIMULATION_BUILDER_H

// GRINS
#include "grins/mesh_builder.h"
#include "grins/visualization_factory.h"
#include "grins/qoi_factory.h"
#include "grins/postprocessing_factory.h"
#include <memory>

namespace GRINS
{
  class SimulationBuilder
  {
  public:

    SimulationBuilder();
    virtual ~SimulationBuilder(){};

    std::shared_ptr<libMesh::UnstructuredMesh> build_mesh
    ( const GetPot& input,
      const libMesh::Parallel::Communicator &comm );

    std::shared_ptr<Visualization> build_vis
    ( const GetPot& input,
      const libMesh::Parallel::Communicator &comm );

    std::shared_ptr<CompositeQoI> build_qoi( const GetPot& input );

    std::shared_ptr<PostProcessedQuantities<libMesh::Real> > build_postprocessing( const GetPot& input );

    void attach_mesh_builder( std::shared_ptr<MeshBuilder> mesh_builder );

    void attach_vis_factory( std::shared_ptr<VisualizationFactory> vis_factory );

    void attach_qoi_factory( std::shared_ptr<QoIFactory> qoi_factory );

    void attach_postprocessing_factory( std::shared_ptr<PostprocessingFactory> postprocessing_factory );

    const MeshBuilder& mesh_builder() const;

  protected:

    std::shared_ptr<MeshBuilder> _mesh_builder;
    std::shared_ptr<VisualizationFactory> _vis_factory;
    std::shared_ptr<QoIFactory> _qoi_factory;
    std::shared_ptr<PostprocessingFactory> _postprocessing_factory;

  }; //class SimulationBuilder
} // namespace GRINS

#endif //GRINS_SIMULATION_BUILDER_H
