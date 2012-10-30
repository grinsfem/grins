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

#ifndef GRINS_SIMULATION_H
#define GRINS_SIMULATION_H

#include "boost/tr1/memory.hpp"

// libMesh
#include "getpot.h"

// GRINS
#include "simulation_builder.h"
#include "visualization.h"
#include "boundary_conditions.h"

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
		SimulationBuilder& sim_builder );

    ~Simulation();
	
    void run();

    void print_sim_info();

    std::tr1::shared_ptr<libMesh::EquationSystems> get_equation_system();	      

    Number get_qoi( unsigned int qoi_index ) const;

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

    std::tr1::shared_ptr<QoIBase> _qoi;

    // Screen display options
    bool _print_mesh_info;
    bool _print_log_info;
    bool _print_equation_system_info;
    bool _print_qoi;

    // Visualization options
    bool _output_vis;
    bool _output_residual;

  };
}
#endif // GRINS_SIMULATION_H
