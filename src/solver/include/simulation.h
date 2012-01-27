//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010,2011 The PECOS Development Team
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

#ifndef SIMULATION_H
#define SIMULATION_H

#include <memory>

// libMesh
#include "getpot.h"

// GRINS
#include "physics_factory.h"
#include "mesh_builder.h"
#include "solver_factory.h"
#include "grins_solver.h"
#include "visualization_factory.h"
#include "visualization.h"

namespace GRINS
{
  class Simulation
  {
  public:
    
    Simulation( const GetPot& input,
		GRINS::PhysicsFactory* physics_factory,
		GRINS::MeshBuilder* mesh_builder,
		GRINS::SolverFactory* solver_factory,
		GRINS::VisualizationFactory* vis_factory );

    ~Simulation();
	
    void run();

    void print_sim_info();

  private:

    std::tr1::shared_ptr<libMesh::Mesh> _mesh;

    std::tr1::shared_ptr<libMesh::EquationSystems> _equation_system;

    std::tr1::shared_ptr<GRINS::Solver> _solver;

    std::tr1::shared_ptr<GRINS::Visualization> _vis;

    // Screen display options
    bool _print_mesh_info;
    bool _print_log_info;
    bool _print_equation_system_info;

    // Visualization options
    bool _output_vis;
    bool _output_residual;

  };
}
#endif // SIMULATION_H
