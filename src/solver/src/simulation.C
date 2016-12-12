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


// This class
#include "grins/simulation.h"

// GRINS
#include "grins/grins_enums.h"
#include "grins/simulation_builder.h"
#include "grins/multiphysics_sys.h"
#include "grins/solver_context.h"
#include "grins/simulation_parsing.h"
#include "grins/strategies_parsing.h"
#include "grins/physics_builder.h"
#include "grins/error_estimator_factory_base.h"
#include "grins/variable_builder.h"

// libMesh
#include "libmesh/dof_map.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/qoi_set.h"
#include "libmesh/sensitivity_data.h"

namespace GRINS
{

  Simulation::Simulation( const GetPot& input,
                          SimulationBuilder& sim_builder,
                          const libMesh::Parallel::Communicator &comm )
    :  _mesh( sim_builder.build_mesh(input, comm) ),
       _equation_system( new libMesh::EquationSystems( *_mesh ) ),
       _solver( sim_builder.build_solver(input) ),
       _system_name( input("screen-options/system_name", "GRINS" ) ),
       _multiphysics_system( &(_equation_system->add_system<MultiphysicsSystem>( _system_name )) ),
       _vis( sim_builder.build_vis(input, comm) ),
       _postprocessing( sim_builder.build_postprocessing(input) ),
       _print_mesh_info( input("screen-options/print_mesh_info", false ) ),
       _print_log_info( input("screen-options/print_log_info", false ) ),
       _print_equation_system_info( input("screen-options/print_equation_system_info", false ) ),
       _print_constraint_info( input("screen-options/print_constraint_info", false ) ),
       _print_scalars( input("screen-options/print_scalars", false ) ),
       _qoi_output( new QoIOutput(input) ),
       _output_vis( input("vis-options/output_vis", false ) ),
       _output_adjoint( input("vis-options/output_adjoint", false ) ),
       _output_residual( input( "vis-options/output_residual", false ) ),
       _output_residual_sensitivities( input( "vis-options/output_residual_sensitivities", false ) ),
       _output_solution_sensitivities( input( "vis-options/output_solution_sensitivities", false ) ),
       _timesteps_per_vis( input("vis-options/timesteps_per_vis", 1 ) ),
       _timesteps_per_perflog( input("screen-options/timesteps_per_perflog", 0 ) ),
       _error_estimator_options(input),
       _error_estimator(), // effectively NULL
       _do_adjoint_solve(false), // Helper function will set final value
       _have_restart(false)
  {
    libmesh_deprecated();

    this->init_multiphysics_system(input);

    this->init_qois(input,sim_builder);

    this->init_params(input,sim_builder);

    this->init_adjoint_solve(input,_output_adjoint);

    // Must be called after setting QoI on the MultiphysicsSystem
    this->build_error_estimator(input);

    if( SimulationParsing::have_restart(input) )
        this->init_restart(input,sim_builder,comm);

    this->check_for_unused_vars(input, false /*warning only*/);

  }

  Simulation::Simulation( const GetPot& input,
                          GetPot& command_line,
                          SimulationBuilder& sim_builder,
                          const libMesh::Parallel::Communicator &comm )
    :  _mesh( sim_builder.build_mesh(input, comm) ),
       _equation_system( new libMesh::EquationSystems( *_mesh ) ),
       _solver( sim_builder.build_solver(input) ),
       _system_name( input("screen-options/system_name", "GRINS" ) ),
       _multiphysics_system( &(_equation_system->add_system<MultiphysicsSystem>( _system_name )) ),
       _vis( sim_builder.build_vis(input, comm) ),
       _postprocessing( sim_builder.build_postprocessing(input) ),
       _print_mesh_info( input("screen-options/print_mesh_info", false ) ),
       _print_log_info( input("screen-options/print_log_info", false ) ),
       _print_equation_system_info( input("screen-options/print_equation_system_info", false ) ),
       _print_constraint_info( input("screen-options/print_constraint_info", false ) ),
       _print_scalars( input("screen-options/print_scalars", false ) ),
       _qoi_output( new QoIOutput(input) ),
       _output_vis( input("vis-options/output_vis", false ) ),
       _output_adjoint( input("vis-options/output_adjoint", false ) ),
       _output_residual( input( "vis-options/output_residual", false ) ),
       _output_residual_sensitivities( input( "vis-options/output_residual_sensitivities", false ) ),
       _output_solution_sensitivities( input( "vis-options/output_solution_sensitivities", false ) ),
       _timesteps_per_vis( input("vis-options/timesteps_per_vis", 1 ) ),
       _timesteps_per_perflog( input("screen-options/timesteps_per_perflog", 0 ) ),
       _error_estimator_options(input),
       _error_estimator(), // effectively NULL
       _do_adjoint_solve(false), // Helper function will set final value
       _have_restart(false)
  {
    this->init_multiphysics_system(input);

    this->init_qois(input,sim_builder);

    this->init_params(input,sim_builder);

    this->init_adjoint_solve(input,_output_adjoint);

    // Must be called after setting QoI on the MultiphysicsSystem
    this->build_error_estimator(input);

    if( SimulationParsing::have_restart(input) )
        this->init_restart(input,sim_builder,comm);

    bool warning_only = command_line.search("--warn-only-unused-var");
    this->check_for_unused_vars(input, warning_only );

  }

  void Simulation::init_multiphysics_system( const GetPot& input )
  {
    // Only print libMesh logging info if the user requests it
    libMesh::perflog.disable_logging();
    if( this->_print_log_info ) libMesh::perflog.enable_logging();

    VariableBuilder::build_variables(input,(*_multiphysics_system));

    PhysicsList physics_list = PhysicsBuilder::build_physics_map(input);

    _multiphysics_system->attach_physics_list( physics_list );

    _multiphysics_system->read_input_options( input );

    _multiphysics_system->register_postprocessing_vars( input, *(_postprocessing) );

    /* Postprocessing needs to be initialized before the solver since that's
       where equation_system gets init'ed */
    _postprocessing->initialize( *_multiphysics_system, *_equation_system );

    _solver->initialize( input, _equation_system, _multiphysics_system );

    // Useful for debugging
    if( input("screen-options/print_dof_constraints", false ) )
      {
        _multiphysics_system->get_dof_map().print_dof_constraints();
      }

    // This *must* be done after equation_system->init in order to get variable indices
    // Set any extra quadrature order the user requested. By default, is 0.
    _multiphysics_system->extra_quadrature_order = StrategiesParsing::extra_quadrature_order(input);
  }

  void Simulation::init_qois( const GetPot& input, SimulationBuilder& sim_builder )
  {
    // If the user actually asks for a QoI, then we add it.
    SharedPtr<CompositeQoI> qois = sim_builder.build_qoi( input );
    if( qois->n_qois() > 0 )
      {
        // This *must* be done after equation_system->init in order to get variable indices
        qois->init(input, *_multiphysics_system );

        /* Note that we are effectively transfering ownership of the qoi pointer because
           it will be cloned in _multiphysics_system and all the calculations are done there. */
        _multiphysics_system->attach_qoi( qois.get() );
      }
    else if (_qoi_output->output_qoi_set())
      {
        std::cout << "Error: print_qoi is specified but\n" <<
          "no QoIs have been specified.\n" << std::endl;
        libmesh_error();
      }

    return;
  }

  void Simulation::init_params( const GetPot& input,
                                SimulationBuilder& /*sim_builder*/ )
  {
    unsigned int n_adjoint_parameters =
      input.vector_variable_size("QoI/adjoint_sensitivity_parameters");

    unsigned int n_forward_parameters =
      input.vector_variable_size("QoI/forward_sensitivity_parameters");

    // If the user actually asks for parameter sensitivities, then we
    // set up the parameter vectors to use.
    if ( n_adjoint_parameters )
      {
        // If we're doing adjoint sensitivities, dq/dp only makes
        // sense if we have q
        CompositeQoI* qoi =
          libMesh::cast_ptr<CompositeQoI*>
            (this->_multiphysics_system->get_qoi());

        if (!qoi)
          {
            std::cout <<
              "Error: adjoint_sensitivity_parameters are specified but\n"
              << "no QoIs have been specified.\n" << std::endl;
            libmesh_error();
          }

        _adjoint_parameters.initialize
          (input, "QoI/adjoint_sensitivity_parameters",
           *this->_multiphysics_system, qoi);
      }

    if ( n_forward_parameters )
      {
        // If we're doing forward sensitivities, du/dp can make
        // sense even with no q defined
        CompositeQoI* qoi =
          dynamic_cast<CompositeQoI*>
            (this->_multiphysics_system->get_qoi());

        // dynamic_cast returns NULL if our QoI isn't a CompositeQoI;
        // i.e. if there were no QoIs that made us bother setting up
        // the CompositeQoI object.  Passing NULL tells
        // ParameterManager not to bother asking for qoi registration
        // of parameters.

        _forward_parameters.initialize
          (input, "QoI/forward_sensitivity_parameters",
           *this->_multiphysics_system, qoi);
      }
  }


  void Simulation::init_restart( const GetPot& input, SimulationBuilder& sim_builder,
                                 const libMesh::Parallel::Communicator &comm )
  {
    this->read_restart( input );

    /* We do this here only if there's a restart file. Otherwise, this was done
       at mesh construction time */
    sim_builder.mesh_builder().do_mesh_refinement_from_input( input, comm, *_mesh );

    /* \todo Any way to tell if the mesh got refined so we don't unnecessarily
       call reinit()? */
    _equation_system->reinit();

    return;
  }

  void Simulation::check_for_unused_vars( const GetPot& input, bool warning_only )
  {
    /* Everything should be set up now, so check if there's any unused variables
       in the input file. If so, then tell the user what they were and error out. */
    std::vector<std::string> unused_vars = input.unidentified_variables();

    bool unused_vars_detected = false;

    // String we use to check if there is a variable in the Materials section
    std::string mat_string = "Materials/";

    // If we do have unused variables, make sure they are in the Materials
    // section since we're allowing those to be present and not used.
    if( !unused_vars.empty() )
      {
        int n_total = unused_vars.size();

        int n_mats = 0;

        for( int v = 0; v < n_total; v++ )
            if( (unused_vars[v]).find(mat_string) != std::string::npos )
              n_mats += 1;

        libmesh_assert_greater_equal( n_total, n_mats );

        if( n_mats < n_total )
          unused_vars_detected = true;
      }

    if( unused_vars_detected )
      {
        libMesh::err << "==========================================================" << std::endl;
        if( warning_only )
          libMesh::err << "Warning: ";
        else
          libMesh::err << "Error: ";

        libMesh::err << "Found unused variables!" << std::endl;

        for( std::vector<std::string>::const_iterator it = unused_vars.begin();
             it != unused_vars.end(); ++it )
          {
            // Don't print out any unused Material variables, since we allow these
            if( (*it).find(mat_string) != std::string::npos )
              continue;

            libMesh::err << *it << std::endl;
          }
        libMesh::err << "==========================================================" << std::endl;

        if( !warning_only )
          libmesh_error();
      }

    return;
  }

  void Simulation::run()
  {
    this->print_sim_info();

    SolverContext context;
    context.system = _multiphysics_system;
    context.equation_system = _equation_system;
    context.vis = _vis;
    context.output_adjoint = _output_adjoint;
    context.timesteps_per_vis = _timesteps_per_vis;
    context.timesteps_per_perflog = _timesteps_per_perflog;
    context.output_vis = _output_vis;
    context.output_residual = _output_residual;
    context.output_residual_sensitivities = _output_residual_sensitivities;
    context.output_solution_sensitivities = _output_solution_sensitivities;
    context.print_scalars = _print_scalars;
    context.print_perflog = _print_log_info;
    context.postprocessing = _postprocessing;
    context.error_estimator = _error_estimator;
    context.qoi_output = _qoi_output;
    context.do_adjoint_solve = _do_adjoint_solve;
    context.have_restart = _have_restart;

    if (_output_residual_sensitivities &&
        !_forward_parameters.parameter_vector.size())
    {
      std::cout <<
        "Error: output_residual_sensitivities is specified but\n" <<
        "no forward sensitivity parameters have been specified.\n" <<
        std::endl;
      libmesh_error();
    }

    if (_output_solution_sensitivities &&
        !_forward_parameters.parameter_vector.size())
    {
      std::cout <<
        "Error: output_solution_sensitivities is specified but\n" <<
        "no forward sensitivity parameters have been specified.\n" <<
        std::endl;
      libmesh_error();
    }

    _solver->solve( context );

    if (_qoi_output->output_qoi_set())
      {
        _multiphysics_system->assemble_qoi();
        const CompositeQoI * my_qoi = libMesh::cast_ptr<const CompositeQoI*>(this->_multiphysics_system->get_qoi());
        _qoi_output->output_qois(*my_qoi, this->_multiphysics_system->comm() );
      }

    if ( _adjoint_parameters.parameter_vector.size() )
      {
        // Default: "calculate sensitivities for all QoIs"
        libMesh::QoISet qois;

        const libMesh::ParameterVector & params =
          _adjoint_parameters.parameter_vector;

        libMesh::SensitivityData sensitivities
          (qois, *this->_multiphysics_system, params);

        _solver->adjoint_qoi_parameter_sensitivity
          (context, qois, params, sensitivities);

        std::cout << "Adjoint sensitivities:" << std::endl;

        for (unsigned int q=0;
             q != this->_multiphysics_system->qoi.size(); ++q)
          {
            for (unsigned int p=0; p != params.size(); ++p)
              {
                std::cout << "dq" << q << "/dp" << p << " = " <<
                        sensitivities[q][p] << std::endl;
              }
          }
      }

    if ( _forward_parameters.parameter_vector.size() )
      {
        // Default: "calculate sensitivities for all QoIs"
        libMesh::QoISet qois;

        const libMesh::ParameterVector & params =
          _forward_parameters.parameter_vector;

        libMesh::SensitivityData sensitivities
          (qois, *this->_multiphysics_system, params);

        _solver->forward_qoi_parameter_sensitivity
          (context, qois, params, sensitivities);

        std::cout << "Forward sensitivities:" << std::endl;

        for (unsigned int q=0;
             q != this->_multiphysics_system->qoi.size(); ++q)
          {
            for (unsigned int p=0; p != params.size(); ++p)
              {
                std::cout << "dq" << q << "/dp" << p << " = " <<
                        sensitivities[q][p] << std::endl;
              }
          }
      }

    return;
  }

  void Simulation::print_sim_info()
  {
    // Print mesh info if the user wants it
    if( this->_print_mesh_info )
      this->_mesh->print_info();

    // Print EquationSystems info if requested
    if( this->_print_equation_system_info )
      this->_equation_system->print_info();

    // Print DofMap constraint info if requested
    if( this->_print_constraint_info )
      this->_multiphysics_system->get_dof_map().print_dof_constraints();

    return;
  }

  SharedPtr<libMesh::EquationSystems> Simulation::get_equation_system()
  {
    return _equation_system;
  }

  libMesh::Number Simulation::get_qoi_value( unsigned int qoi_index ) const
  {
    const CompositeQoI* qoi = libMesh::cast_ptr<const CompositeQoI*>(this->_multiphysics_system->get_qoi());
    return qoi->get_qoi_value(qoi_index);
  }

  void Simulation::read_restart( const GetPot& input )
  {
    const std::string restart_file = SimulationParsing::restart_file(input);

    // Most of this was pulled from FIN-S
    if (restart_file != "none")
      {
        std::cout << " ====== Restarting from " << restart_file << std::endl;

        // Must have correct file type to restart
        if (restart_file.rfind(".xdr") < restart_file.size())
          _equation_system->read(restart_file,GRINSEnums::DECODE,
                                 //libMesh::EquationSystems::READ_HEADER |  // Allow for thermochemistry upgrades
                                 libMesh::EquationSystems::READ_DATA |
                                 libMesh::EquationSystems::READ_ADDITIONAL_DATA);

        else if  (restart_file.rfind(".xda") < restart_file.size())
          _equation_system->read(restart_file,GRINSEnums::READ,
                                 //libMesh::EquationSystems::READ_HEADER |  // Allow for thermochemistry upgrades
                                 libMesh::EquationSystems::READ_DATA |
                                 libMesh::EquationSystems::READ_ADDITIONAL_DATA);

        else
          {
            std::cerr << "Error: Restart filename must have .xdr or .xda extension!" << std::endl;
            libmesh_error();
          }

        this->_have_restart = true;

        const std::string system_name = input("screen-options/system_name", "GRINS" );

        MultiphysicsSystem& system =
          _equation_system->get_system<MultiphysicsSystem>(system_name);

        // Update the old data
        system.update();
      }

    return;
  }

  void Simulation::init_adjoint_solve( const GetPot& input, bool output_adjoint )
  {
    // Check if we're doing an adjoint solve
    _do_adjoint_solve = this->check_for_adjoint_solve(input);

    const libMesh::DifferentiableQoI* raw_qoi = _multiphysics_system->get_qoi();
    const CompositeQoI* qoi = dynamic_cast<const CompositeQoI*>( raw_qoi );

    // If we are trying to do an adjoint solve without a QoI, that's an error
    // If there are no QoIs, the CompositeQoI never gets built and qoi will be NULL
    if( _do_adjoint_solve && !qoi )
      {
        libmesh_error_msg("Error: Adjoint solve requested, but no QoIs detected.");
      }

    /* If we are not doing an adjoint solve, but adjoint output was requested:
       that's an error. */
    if( !_do_adjoint_solve && output_adjoint )
      {
        libmesh_error_msg("Error: Adjoint output requested, but no adjoint solve requested.");
      }
  }

  bool Simulation::check_for_adjoint_solve( const GetPot& input ) const
  {
    std::string error_estimator = _error_estimator_options.estimator_type();

    bool do_adjoint_solve = false;

    // First check if the error estimator requires an adjoint solve
    if( _error_estimator_options.estimator_requires_adjoint() )
      do_adjoint_solve = true;

    // Next check if parameter sensitivies require an adjoint solve
    if ( _adjoint_parameters.parameter_vector.size() )
      do_adjoint_solve = true;

    // Now check if the user requested to do an adjoint solve regardless
    /*! \todo This option name WILL change at some point. */
    if( StrategiesParsing::do_adjoint_solve(input) )
      do_adjoint_solve = true;

    return do_adjoint_solve;
  }

  void Simulation::build_error_estimator(const GetPot& input)
  {
    std::string estimator_type = _error_estimator_options.estimator_type();
    if( estimator_type != std::string("none") )
      {
        ErrorEstimatorFactoryBase::set_getpot(input);
        ErrorEstimatorFactoryBase::set_system(*_multiphysics_system);
        ErrorEstimatorFactoryBase::set_estimator_options(_error_estimator_options);
        _error_estimator.reset( (ErrorEstimatorFactoryBase::build(estimator_type)).release() );
      }
  }

#ifdef GRINS_USE_GRVY_TIMERS
  void Simulation::attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer )
  {
    this->_multiphysics_system->attach_grvy_timer( grvy_timer );
    return;
  }
#endif

} // namespace GRINS
