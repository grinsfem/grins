//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "grins/simulation_builder.h"

namespace GRINS
{
  SimulationBuilder::SimulationBuilder()
    : _physics_factory( new PhysicsFactory ),
      _mesh_builder( new MeshBuilder ),
      _solver_factory( new SolverFactory ),
      _vis_factory( new VisualizationFactory ),
      _bc_factory( new BoundaryConditionsFactory ),
      _qoi_factory( new QoIFactory ),
      _postprocessing_factory( new PostprocessingFactory )
  {
    return;
  }

  SimulationBuilder::~SimulationBuilder()
  {
    return;
  }
  
  void SimulationBuilder::attach_physics_factory( std::tr1::shared_ptr<PhysicsFactory> physics_factory )
  {
    this->_physics_factory = physics_factory;
    return;
  }
  
  void SimulationBuilder::attach_solver_factory( std::tr1::shared_ptr<SolverFactory> solver_factory )
  {
    this->_solver_factory = solver_factory;
    return;
  }

  void SimulationBuilder::attach_mesh_builder( std::tr1::shared_ptr<MeshBuilder> mesh_builder )
  {
    this->_mesh_builder = mesh_builder;
    return;
  }
    
  void SimulationBuilder::attach_vis_factory( std::tr1::shared_ptr<VisualizationFactory> vis_factory )
  {
    this->_vis_factory = vis_factory;
    return;
  }

  void SimulationBuilder::attach_bc_factory( std::tr1::shared_ptr<BoundaryConditionsFactory> bc_factory )
  {
    this->_bc_factory = bc_factory;
    return;
  }

  void SimulationBuilder::attach_qoi_factory( std::tr1::shared_ptr<QoIFactory> qoi_factory )
  {
    this->_qoi_factory = qoi_factory;
  }

  void SimulationBuilder::attach_postprocessing_factory( std::tr1::shared_ptr<PostprocessingFactory> postprocessing_factory )
  {
    this->_postprocessing_factory = postprocessing_factory; 
  }

  std::tr1::shared_ptr<libMesh::Mesh> SimulationBuilder::build_mesh( const GetPot& input )
  {
    return (this->_mesh_builder)->build(input);
  }

  GRINS::PhysicsList SimulationBuilder::build_physics( const GetPot& input )
  {
    return (this->_physics_factory)->build(input);
  }
  
  std::tr1::shared_ptr<GRINS::Solver> SimulationBuilder::build_solver( const GetPot& input )
  {
    return (this->_solver_factory)->build(input);
  }

    std::tr1::shared_ptr<GRINS::Visualization> SimulationBuilder::build_vis( const GetPot& input )
    {
      return (this->_vis_factory)->build(input);
    }

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > SimulationBuilder::build_dirichlet_bcs()
  {
    return (this->_bc_factory)->build_dirichlet();
  }

  std::map< GRINS::PhysicsName, GRINS::NBCContainer > SimulationBuilder::build_neumann_bcs( libMesh::EquationSystems& equation_system )
  {
    return (this->_bc_factory)->build_neumann(equation_system);
  }
  
  std::tr1::shared_ptr<QoIBase> SimulationBuilder::build_qoi( const GetPot& input )
  {
    return (this->_qoi_factory)->build(input);
  }

  std::tr1::shared_ptr<PostProcessedQuantities<Real> > SimulationBuilder::build_postprocessing( const GetPot& input )
  {
    return (this->_postprocessing_factory)->build(input);
  }

} //namespace GRINS
