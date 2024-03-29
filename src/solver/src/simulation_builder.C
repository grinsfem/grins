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
#include "grins/simulation_builder.h"

// libMesh
#include "libmesh/error_estimator.h"
#include "libmesh/adjoint_refinement_estimator.h"

namespace GRINS
{
  SimulationBuilder::SimulationBuilder()
    : _mesh_builder( std::make_shared<MeshBuilder>() ),
      _vis_factory( std::make_shared<VisualizationFactory>() ),
      _qoi_factory( std::make_shared<QoIFactory>() ),
      _postprocessing_factory( std::make_shared<PostprocessingFactory>() )
  {}

  void SimulationBuilder::attach_mesh_builder( std::shared_ptr<MeshBuilder> mesh_builder )
  {
    this->_mesh_builder = mesh_builder;
  }

  void SimulationBuilder::attach_vis_factory( std::shared_ptr<VisualizationFactory> vis_factory )
  {
    this->_vis_factory = vis_factory;
  }

  void SimulationBuilder::attach_qoi_factory( std::shared_ptr<QoIFactory> qoi_factory )
  {
    this->_qoi_factory = qoi_factory;
  }

  void SimulationBuilder::attach_postprocessing_factory( std::shared_ptr<PostprocessingFactory> postprocessing_factory )
  {
    this->_postprocessing_factory = postprocessing_factory;
  }

  std::shared_ptr<libMesh::UnstructuredMesh> SimulationBuilder::build_mesh
  ( const GetPot& input,
    const libMesh::Parallel::Communicator &comm)
  {
    return (this->_mesh_builder)->build(input, comm);
  }

  std::shared_ptr<Visualization> SimulationBuilder::build_vis
  ( const GetPot& input,
    const libMesh::Parallel::Communicator &comm)
  {
    return (this->_vis_factory)->build(input, comm);
  }

  std::shared_ptr<CompositeQoI> SimulationBuilder::build_qoi( const GetPot& input )
  {
    return (this->_qoi_factory)->build(input);
  }

  std::shared_ptr<PostProcessedQuantities<libMesh::Real> >
  SimulationBuilder::build_postprocessing( const GetPot& input )
  {
    return (this->_postprocessing_factory)->build(input);
  }

  const MeshBuilder& SimulationBuilder::mesh_builder() const
  {
    return *_mesh_builder;
  }

} //namespace GRINS
