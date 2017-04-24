//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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
#include "grins/shared_ptr.h"

// GRINS
#include "grins_config.h"
#include "grins/shared_ptr.h"
#include "grins/grins_solver.h"
#include "grins/qoi_base.h"
#include "grins/visualization.h"
#include "grins/physics_naming.h"
#include "grins/parameter_manager.h"
#include "grins/postprocessed_quantities.h"
#include "grins/error_estimator_options.h"
#include "grins/qoi_output.h"

// libMesh
#include "libmesh/error_estimator.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"

// GRVY
#ifdef GRINS_HAVE_GRVY
#include "libmesh/ignore_warnings.h" // avoid auto_ptr deprecated warnings
#include "grvy.h"
#include "libmesh/restore_warnings.h"
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

    Simulation( const GetPot& input,
                GetPot& command_line, /* Has to be non-const for search() */
		SimulationBuilder& sim_builder,
                const libMesh::Parallel::Communicator &comm
                LIBMESH_CAN_DEFAULT_TO_COMMWORLD );

    virtual ~Simulation(){};

    void run();

    void print_sim_info();

    SharedPtr<libMesh::EquationSystems> get_equation_system();
    MultiphysicsSystem* get_multiphysics_system();

    const MultiphysicsSystem* get_multiphysics_system() const;

    libMesh::Number get_qoi_value( unsigned int qoi_index ) const;

    const std::string& get_multiphysics_system_name() const;

#ifdef GRINS_USE_GRVY_TIMERS
    void attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer );
#endif

  protected:

    void read_restart( const GetPot& input );

    //! Helper function
    void init_multiphysics_system( const GetPot& input );

    //! Helper function
    void init_qois( const GetPot& input, SimulationBuilder& sim_builder );

    //! Helper function
    void init_params( const GetPot& input, SimulationBuilder& sim_builder );

    //! Helper function
    void init_restart( const GetPot& input, SimulationBuilder& sim_builder,
                       const libMesh::Parallel::Communicator &comm );

    //! Helper function
    void check_for_unused_vars( const GetPot& input, bool warning_only );

    //! Helper function
    bool check_for_adjoint_solve( const GetPot& input ) const;

    //! Helper function
    void init_adjoint_solve( const GetPot& input, bool output_adjoint );

    void build_error_estimator(const GetPot& input);

    SharedPtr<libMesh::UnstructuredMesh> _mesh;

    SharedPtr<libMesh::EquationSystems> _equation_system;

    SharedPtr<GRINS::Solver> _solver;

    //! GRINS::Multiphysics system name
    std::string _system_name;
    
    // This needs to be a standard pointer, as _equation_system will own and destroy the object.
    GRINS::MultiphysicsSystem* _multiphysics_system;

    SharedPtr<GRINS::Visualization> _vis;

    SharedPtr<PostProcessedQuantities<libMesh::Real> > _postprocessing;

    // Screen display options
    bool _print_mesh_info;
    bool _print_log_info;
    bool _print_equation_system_info;
    bool _print_constraint_info;
    bool _print_perflog;
    bool _print_scalars;

    // QoI output options and functionality
    SharedPtr<QoIOutput> _qoi_output;

    // Visualization options
    bool _output_vis;
    bool _output_adjoint;
    bool _output_residual;
    bool _output_residual_sensitivities;
    bool _output_solution_sensitivities;

    unsigned int _timesteps_per_vis;
    unsigned int _timesteps_per_perflog;

    ErrorEstimatorOptions _error_estimator_options;
    SharedPtr<libMesh::ErrorEstimator> _error_estimator;

    ParameterManager _adjoint_parameters;

    ParameterManager _forward_parameters;

    // Cache whether or not we do an adjoint solve
    bool _do_adjoint_solve;

    bool _have_restart;

  private:

    Simulation();

  };

  inline
  const MultiphysicsSystem*
  Simulation::get_multiphysics_system() const
  {
    return this->_multiphysics_system;
  }

  inline
  MultiphysicsSystem*
  Simulation::get_multiphysics_system()
  {
    return this->_multiphysics_system;
  }

  inline
  const std::string& Simulation::get_multiphysics_system_name() const
  {
    return this->_system_name;
  }
}
#endif // GRINS_SIMULATION_H
