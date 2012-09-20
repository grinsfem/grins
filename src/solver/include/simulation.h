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

#ifndef SIMULATION_H
#define SIMULATION_H

#include "boost/tr1/memory.hpp"

// libMesh
#include "getpot.h"

// GRINS
#include "physics_factory.h"
#include "mesh_builder.h"
#include "solver_factory.h"
#include "grins_solver.h"
#include "visualization_factory.h"
#include "visualization.h"
#include "boundary_conditions.h"
#include "bc_factory.h"

// GRVY
#ifdef HAVE_GRVY
#include "grvy.h"
#endif

namespace GRINS
{
  class Simulation
  {
  public:
    
    Simulation( const GetPot& input,
		GRINS::PhysicsFactory* physics_factory,
		GRINS::MeshBuilder* mesh_builder,
		GRINS::SolverFactory* solver_factory,
		GRINS::VisualizationFactory* vis_factory,
		GRINS::BoundaryConditionsFactory* bc_factory = NULL );

    ~Simulation();
	
    void run();

    void print_sim_info();

    std::tr1::shared_ptr<libMesh::EquationSystems> get_equation_system();	      

#ifdef USE_GRVY_TIMERS
    void attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer );
#endif

  private:
    
    void check_for_restart( const GetPot& input );

    void attach_neumann_bc_funcs( std::map< GRINS::PhysicsName, GRINS::NBCContainer > neumann_bcs,
				  GRINS::MultiphysicsSystem* system );
    
    void attach_dirichlet_bc_funcs( std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > dbc_map,
				    GRINS::MultiphysicsSystem* system );

    std::tr1::shared_ptr<libMesh::Mesh> _mesh;

    std::tr1::shared_ptr<libMesh::EquationSystems> _equation_system;

    std::tr1::shared_ptr<GRINS::Solver> _solver;

    //! GRINS::Multiphysics system name
    std::string _system_name;
    
    GRINS::MultiphysicsSystem* _multiphysics_system;

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
