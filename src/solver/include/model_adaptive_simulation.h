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
// $Id: mesh_adaptive_simulation.h tvanopstal
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_MODEL_ADAPTIVE_SIMULATION_H
#define GRINS_MODEL_ADAPTIVE_SIMULATION_H

#include "boost/tr1/memory.hpp"

// libMesh
#include "adjoint_refinement_estimator.h"
#include "getpot.h"

// GRINS
#include "boundary_conditions.h"
#include "simulation_builder.h"
#include "visualization.h"

// GRVY
#ifdef HAVE_GRVY
#include "grvy.h"
#endif

namespace GRINS
{
  class ModelAdaptiveSimulation
  {
  public:
    
    ModelAdaptiveSimulation(
        const GetPot& forward_input,
        const GetPot& adjoint_input,
        const GetPot& residual_input,
        SimulationBuilder& forward_sim_builder,
        SimulationBuilder& adjoint_sim_builder,
        SimulationBuilder& residual_sim_builder );

    ~ModelAdaptiveSimulation();
	
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

    // Objects for forward problem
    std::tr1::shared_ptr<libMesh::Mesh> _forward_mesh;
    std::tr1::shared_ptr<libMesh::EquationSystems> _forward_equation_system;
    std::tr1::shared_ptr<GRINS::Solver> _forward_solver;
    std::string _forward_system_name;
    // This needs to be a standard pointer, as _equation_system will own and destroy the object.
    GRINS::MultiphysicsSystem* _forward_multiphysics_system;
    std::tr1::shared_ptr<GRINS::Visualization> _forward_vis;
    std::tr1::shared_ptr<QoIBase> _forward_qoi;

    // Objects for adjoint problem
    std::tr1::shared_ptr<libMesh::Mesh> _adjoint_mesh;
    std::tr1::shared_ptr<libMesh::EquationSystems> _adjoint_equation_system;
    std::tr1::shared_ptr<GRINS::Solver> _adjoint_solver;
    std::string _adjoint_system_name;
    // This needs to be a standard pointer, as _equation_system will own and destroy the object.
    GRINS::MultiphysicsSystem* _adjoint_multiphysics_system;
    std::tr1::shared_ptr<GRINS::Visualization> _adjoint_vis;
    std::tr1::shared_ptr<QoIBase> _adjoint_qoi;

    // Objects for residual
    std::tr1::shared_ptr<libMesh::Mesh> _residual_mesh;
    std::tr1::shared_ptr<libMesh::EquationSystems> _residual_equation_system;
    std::tr1::shared_ptr<GRINS::Solver> _residual_solver;
    std::string _residual_system_name;
    // This needs to be a standard pointer, as _equation_system will own and destroy the object.
    GRINS::MultiphysicsSystem* _residual_multiphysics_system;
    std::tr1::shared_ptr<GRINS::Visualization> _residual_vis;

    // Adaptivity parameters
    unsigned int _max_r_steps;
    bool _output_adjoint_sol;
    double _absolute_global_tolerance;
    bool _plot_cell_errors;
    std::string _error_plot_prefix;

    // Screen display options
    bool _print_mesh_info;
    bool _print_log_info;
    bool _print_equation_system_info;
    bool _print_qoi;

    // Visualization options
    bool _output_vis;
    bool _output_residual;

    std::tr1::shared_ptr<libMesh::AdjointRefinementEstimator> _adjoint_refinement_estimator;

  };
}
#endif // GRINS_MODEL_ADAPTIVE_SIMULATION_H
