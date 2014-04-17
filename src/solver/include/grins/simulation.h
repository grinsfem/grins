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


#ifndef GRINS_SIMULATION_H
#define GRINS_SIMULATION_H

// C++
#include "boost/tr1/memory.hpp"

// GRINS
#include "grins_config.h"
#include "grins/grins_solver.h"
#include "grins/qoi_base.h"
#include "grins/visualization.h"
#include "grins/grins_physics_names.h"
#include "grins/nbc_container.h"
#include "grins/dbc_container.h"
#include "grins/postprocessed_quantities.h"

// libMesh
#include "libmesh/error_estimator.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"

// GRVY
#ifdef GRINS_HAVE_GRVY
#include "grvy.h"
#endif

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  // Forward declarations
  class SimulationBuilder;
  class MultiphysicsSystem;

  class Simulation
  {
  public:
    
    Simulation( const GetPot& input,
		SimulationBuilder& sim_builder,
                const libMesh::Parallel::Communicator &comm 
                LIBMESH_CAN_DEFAULT_TO_COMMWORLD );

    virtual ~Simulation();
	
    void run();

    void print_sim_info();

    std::tr1::shared_ptr<libMesh::EquationSystems> get_equation_system();	      

    libMesh::Number get_qoi_value( unsigned int qoi_index ) const;

#ifdef GRINS_USE_GRVY_TIMERS
    void attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer );
#endif

  protected:
    
    void check_for_restart( const GetPot& input );

    void attach_neumann_bc_funcs( std::map< GRINS::PhysicsName, GRINS::NBCContainer > neumann_bcs,
				  GRINS::MultiphysicsSystem* system );
    
    void attach_dirichlet_bc_funcs( std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > dbc_map,
				    GRINS::MultiphysicsSystem* system );

    std::tr1::shared_ptr<libMesh::UnstructuredMesh> _mesh;

    std::tr1::shared_ptr<libMesh::EquationSystems> _equation_system;

    std::tr1::shared_ptr<GRINS::Solver> _solver;

    //! GRINS::Multiphysics system name
    std::string _system_name;
    
    // This needs to be a standard pointer, as _equation_system will own and destroy the object.
    GRINS::MultiphysicsSystem* _multiphysics_system;

    std::tr1::shared_ptr<GRINS::Visualization> _vis;

    std::tr1::shared_ptr<PostProcessedQuantities<Real> > _postprocessing;

    // Screen display options
    bool _print_mesh_info;
    bool _print_log_info;
    bool _print_equation_system_info;
    bool _print_qoi;

    // Visualization options
    bool _output_vis;
    bool _output_residual;

    unsigned int _timesteps_per_vis;

    std::tr1::shared_ptr<libMesh::ErrorEstimator> _error_estimator;

  private:

    Simulation();

  };
}
#endif // GRINS_SIMULATION_H
