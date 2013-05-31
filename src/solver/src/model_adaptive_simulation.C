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
// $Id: model_adaptive_simulation.C tvanopstal
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "model_error_estimator.h"
#include "model_adaptive_simulation.h"
#include "unistd.h"

GRINS::ModelAdaptiveSimulation::ModelAdaptiveSimulation(
    const GetPot& forward_input,
    const GetPot& adjoint_input,
    const GetPot& residual_input,
	  SimulationBuilder& forward_sim_builder,
	  SimulationBuilder& adjoint_sim_builder,
	  SimulationBuilder& residual_sim_builder )
  : _forward_mesh( forward_sim_builder.build_mesh( forward_input ) ),
    _adjoint_mesh( adjoint_sim_builder.build_mesh( adjoint_input ) ),
    _residual_mesh( residual_sim_builder.build_mesh( residual_input ) ),
    
    _forward_equation_system( new libMesh::EquationSystems( *_forward_mesh ) ),
    _adjoint_equation_system( new libMesh::EquationSystems( *_adjoint_mesh ) ),
    _residual_equation_system( new libMesh::EquationSystems( *_residual_mesh ) ),

    _forward_solver( forward_sim_builder.build_solver( forward_input ) ),
    _adjoint_solver( adjoint_sim_builder.build_solver( adjoint_input ) ),
    _residual_solver( residual_sim_builder.build_solver( residual_input ) ),

    // We won't let user set these for fear of seg faults (systems deleted by system name, names must be unique)
    // _forward_system_name( forward_input("screen-options/system_name", "forward_system" ) ),
    // _adjoint_system_name( adjoint_input("screen-options/system_name", "adjoint_system" ) ),
    // _residual_system_name( residual_input("screen-options/system_name", "residual_system" ) ),

    // System name hard-coded in multiple places in this file!
    // GRINS is default for simulation.C, so we can restart from those simulations with relative ease.
    _forward_multiphysics_system( &(_forward_equation_system->add_system<GRINS::MultiphysicsSystem>( "GRINS" )) ),
    _adjoint_multiphysics_system( &(_adjoint_equation_system->add_system<GRINS::MultiphysicsSystem>( "adjoint_system" )) ),
    _residual_multiphysics_system( &(_residual_equation_system->add_system<GRINS::MultiphysicsSystem>( "residual_system" )) ),

    _forward_vis( forward_sim_builder.build_vis( forward_input ) ),
    _adjoint_vis( adjoint_sim_builder.build_vis( adjoint_input ) ),
    _residual_vis( residual_sim_builder.build_vis( residual_input ) ),

    _forward_qoi( forward_sim_builder.build_qoi( adjoint_input ) ),
    _adjoint_qoi( adjoint_sim_builder.build_qoi( adjoint_input ) ),
    _adjoint_refinement_estimator( adjoint_sim_builder.build_adjoint_refinement_estimator( adjoint_input, _adjoint_qoi ) ),
    // Most output options stashed in forward_input
    _print_mesh_info( forward_input("screen-options/print_mesh_info", false ) ),
    _print_log_info( forward_input("screen-options/print_log_info", false ) ),
    _print_equation_system_info( forward_input("screen-options/print_equation_system_info", false ) ),
    _print_qoi( adjoint_input("screen-options/print_qoi", false ) ),
    _output_vis( forward_input("vis-options/output_vis", false ) ),
    _output_residual( forward_input( "vis-options/output_residual", false ) )
{
  // Only print libMesh logging info if the user requests it
  libMesh::perflog.disable_logging();
  if( this->_print_log_info ) libMesh::perflog.enable_logging();

  // Read input options
  this->_max_r_steps = forward_input( "Adaptivity/max_r_steps", 1 );
  this->_output_adjoint_sol = forward_input( "Adaptivity/output_adjoint_sol", false );
  this-> _absolute_global_tolerance = forward_input( "Adaptivity/absolute_global_tolerance", -1. );
  this->_plot_cell_errors = forward_input( "Adaptivity/plot_cell_errors", false );
  this->_error_plot_prefix = forward_input( "Adaptivity/error_plot_prefix", "cell_error" );

  //-------------------------------------
  // Simulation setup for forward problem
  //-------------------------------------
  out << "\nModel adaptive simulation: forward problem" << std::endl;
  GRINS::PhysicsList forward_physics_list = forward_sim_builder.build_physics( forward_input );
  _forward_multiphysics_system->attach_physics_list( forward_physics_list );
  _forward_multiphysics_system->read_input_options( forward_input );

  // This *must* be done before equation_system->init
  this->attach_dirichlet_bc_funcs( forward_sim_builder.build_dirichlet_bcs(), &*_forward_multiphysics_system );

  // Includes _forward_equation_system->init
  _forward_solver->initialize( forward_input, _forward_equation_system, &*_forward_multiphysics_system );

  // This *must* be done after equation_system->init in order to get variable indices
  this->attach_neumann_bc_funcs( forward_sim_builder.build_neumann_bcs( *_forward_equation_system ), &*_forward_multiphysics_system );

  // Assert that we have at least one QoI, and for the moment not more than one either.
  if( _forward_qoi.use_count() == 1 )
  {
    // This *must* be done after equation_system->init in order to get variable indices
    this->_forward_qoi->init( forward_input, *_forward_multiphysics_system );
    
    _forward_multiphysics_system->attach_qoi( &(*(this->_forward_qoi)) );
  }
  else
  {
    out << "Error: for now, only single-valued QoI sets are accepted." << std::endl;
    libmesh_error();
  }

  //-------------------------------------------
  // Largely identical setup of adjoint problem
  //-------------------------------------------
  out << "\nModel adaptive simulation: adjoint problem" << std::endl;
  GRINS::PhysicsList adjoint_physics_list = adjoint_sim_builder.build_physics( adjoint_input );
  _adjoint_multiphysics_system->attach_physics_list( adjoint_physics_list );
  _adjoint_multiphysics_system->read_input_options( adjoint_input );

  // This *must* be done before equation_system->init
  this->attach_dirichlet_bc_funcs( adjoint_sim_builder.build_dirichlet_bcs(), &*_adjoint_multiphysics_system );

  _adjoint_solver->initialize( adjoint_input, _adjoint_equation_system, &*_adjoint_multiphysics_system );

  // This *must* be done after equation_system->init in order to get variable indices
  this->attach_neumann_bc_funcs( adjoint_sim_builder.build_neumann_bcs( *_adjoint_equation_system ), &*_adjoint_multiphysics_system );

  // Assert that we have at least one QoI, and for the moment not more than one either.
  if( _adjoint_qoi.use_count() == 1 )
  {
    // This *must* be done after equation_system->init in order to get variable indices
    this->_adjoint_qoi->init( adjoint_input, *_adjoint_multiphysics_system );
    
    _adjoint_multiphysics_system->attach_qoi( &(*(this->_adjoint_qoi)) );
  }
  else
  {
    out << "Error: for now, only single-valued QoI sets are accepted." << std::endl;
    libmesh_error();
  }

  //----------------------------------------------------------------------------
  // Residual and DWR refinement (largely identical to setup of forward problem)
  //----------------------------------------------------------------------------
  out << "\nModel adaptive simulation: residual" << std::endl;
  GRINS::PhysicsList residual_physics_list = residual_sim_builder.build_physics( residual_input );
  _residual_multiphysics_system->attach_physics_list( residual_physics_list );
  _residual_multiphysics_system->read_input_options( residual_input );

  // This *must* be done before equation_system->init
  this->attach_dirichlet_bc_funcs( residual_sim_builder.build_dirichlet_bcs(), &*_residual_multiphysics_system );

  _residual_solver->initialize( residual_input, _residual_equation_system, &*_residual_multiphysics_system );

  // This *must* be done after equation_system->init in order to get variable indices
  this->attach_neumann_bc_funcs( residual_sim_builder.build_neumann_bcs( *_residual_equation_system ), &*_residual_multiphysics_system );

  this->check_for_restart( forward_input );

  return;
}

GRINS::ModelAdaptiveSimulation::~ModelAdaptiveSimulation()
{
  return;
}

void GRINS::ModelAdaptiveSimulation::run()
{
  this->print_sim_info();

  // Empty objects
  // std::tr1::shared_ptr<GRINS::QoIBase> empty_qoi = std::tr1::shared_ptr<GRINS::QoIBase>();
  std::tr1::shared_ptr<libMesh::ErrorEstimator> empty_error_estimator =
      std::tr1::shared_ptr<libMesh::ErrorEstimator>();

  // Model adaptive loop
  for( unsigned int r_step = 0;
       r_step != this->_max_r_steps;
       r_step++ )
  {
    //------------------------------------------------------------
    // Solve forward problem (giving empty QoI and ErrorEstimator)
    //------------------------------------------------------------

    out << "\nForward solve...\n";
    _forward_solver->solve( &*_forward_multiphysics_system, _forward_equation_system,
        _forward_qoi, _forward_vis, false, false, empty_error_estimator );

    if( _output_vis )
    {
      out << "@TO: output primal (from primal system)\n";
      // Write out primal solution
      _forward_vis->output( _forward_equation_system );
    }

    if( this->_print_qoi )
    {
      libMesh::QoISet qoi_set( *_forward_multiphysics_system );
      _forward_multiphysics_system->assemble_qoi( qoi_set ) ;
      const QoIBase* qoi = libmesh_cast_ptr<const QoIBase*>(this->_forward_multiphysics_system->get_qoi());
      qoi->output_qoi( std::cout );
    }

    _forward_multiphysics_system->assembly( true, false );
    _forward_multiphysics_system->rhs->close();
    out << "1 @TO: residual (from forward system), with linfty norm "
        << _forward_multiphysics_system->rhs->linfty_norm() << " \n";

    /*
    //----------------------------------------------------
    // Solve adjoint problem (in same approximation space)
    //----------------------------------------------------

    // Substitute primal solution from forward problem into adjoint system
    // For now this does not involve a projection step
    NumericVector<Number>& primal_of_adj = *_adjoint_multiphysics_system->solution;
    NumericVector<Number>& primal = *_forward_multiphysics_system->solution;
    primal_of_adj.swap( primal );

    // Solve
    out << "\nAdjoint solve...\n";
    _adjoint_multiphysics_system->adjoint_solve();

    // Swap back primal solution into primal system
    primal_of_adj.swap( primal ); // Swap back for use in residual system
  
    // Output
    NumericVector<Number>& dual = (_adjoint_multiphysics_system->get_adjoint_solution(0));
    if( _output_vis )
    {
      out << "@TO: output adjoint (from adjoint system)\n";
      // Swap primal and dual to write out dual solution
	    primal_of_adj.swap( dual );	    
      _adjoint_vis->output( _adjoint_equation_system );
	    primal_of_adj.swap( dual );	    
    }
    */

    //--------------------------------
    // Evaluate dual weighted residual
    //--------------------------------

    //_residual_multiphysics_system->reinit();
    out << "\nError estimation and adaptivity...\n";

    /* TODO PUT ME BACK!!!
    // Set forward solution
    NumericVector<Number> &primal_of_res = *_residual_multiphysics_system->solution;
    primal_of_res.swap( primal );

    // Set adjoint solution
    _residual_multiphysics_system->add_adjoint_solution( 0 );
    NumericVector<Number> &dual_of_res = _residual_multiphysics_system->get_adjoint_solution( 0 );
    dual_of_res.swap( dual );
    _forward_multiphysics_system->add_adjoint_solution( 0 );
    */

    // My trash ----------------------------------------------------------------------
    MeshBase::element_iterator iter_begin = _forward_mesh->active_local_elements_begin();
    MeshBase::element_iterator iter_end = _forward_mesh->active_local_elements_end();
    for( libMesh::MeshBase::element_iterator elem_it = iter_begin;
         elem_it != iter_end;
         ++elem_it )
    {
      const Elem *elem = const_cast<Elem *>( *elem_it );
      unsigned int e_id = (*elem_it)->id();
      std::vector<unsigned int> dof_indices;
      _forward_multiphysics_system->get_dof_map().dof_indices( elem, dof_indices );
      unsigned int n_dofs = dof_indices.size();
      // Print out most everything
      out << "(b) elem # " << e_id << "\n";
      std::vector<Number> tmpvec;
      if( true )
      {
        _forward_multiphysics_system->solution->get( dof_indices, tmpvec );
      }
      else
      {
        _forward_multiphysics_system->rhs->get( dof_indices, tmpvec );
      }
      for( unsigned int j = 0; j != tmpvec.size(); ++j )
      {
        out << tmpvec[j] << "\n";
      }
      out << "\n";
    }
    // -------------------------------------------------------------------------------

    // Initialize DWR model error estimator object
    ErrorVector error;
    // ---------------------------------------------------------------------------------------
    // for( ConstElemRange::const_iterator elem_it = _forward_mesh.active_local_elements_begin();
    //      elem_it != _forward_mesh.active_local_elements_end();
    //      ++elem_it )
    // {
    //   const Elem *elem = const_cast<Elem *>( *elem_it );
    //   // Element dof indices and size
    //   std::vector<unsigned int> dof_indices;
    //   _sys->get_dof_map().dof_indices( elem, dof_indices );
    //   // Reference to created element residual
    //   DenseVector<Number> &elem_residual = _fem_context.elem_residual;
    //   unsigned int n_dofs = dof_indices.size();
    //   for( unsigned int i=0; i != n_dofs; ++i )
    //   {
    //     elem_residual(i) = 
    //   }
    //   unsigned int e_id = (*elem_it)->id();
    // ---------------------------------------------------------------------------------------

    // TODO PUT ME BACK!!!
    // ModelErrorEstimator model_err_estimator(
    //     _residual_multiphysics_system,
    //     &_residual_multiphysics_system->get_mesh(),
    //     &error );
    ModelErrorEstimator model_err_estimator(
        _forward_multiphysics_system,
        &*_forward_mesh,
        &error );
    
    // Compute error based on DWR
    model_err_estimator.assembly();

    // Plot error vector
    if( this->_plot_cell_errors )
        error.plot_error( this->_error_plot_prefix+".exo", *(this->_residual_mesh) );

    if( _output_residual )
    {
      // Write out residual at dual solution
      _residual_vis->output_residual( _residual_equation_system, &*_residual_multiphysics_system );
    }

    // After error estimation, solutions can be swapped back
    // TODO PUT ME BACK!!!
    // primal_of_res.swap( primal );
    // dual_of_res.swap( dual );

    /*
    // Flag and ``refine''
    ModelRefinement model_refinement;
	  if( this->_absolute_global_tolerance >= 0. )
	  {
	    model_refinement.flag_elements_by_error_tolerance( error );
	  }
	  else
	  {
	    std::cout << "NotImplementedError: Unknown model refinement strategy."<< std::endl;
	  }
	  model_refinement.relabel_elements(); // implement this
    */

    // Reinit equation_system (necessary in case of mesh adaptive solvers)
	  //_forward_equation_system->reinit();
	  //_adjoint_equation_system->reinit();

    // Output summary
    out << "At model-adaptive step: " << r_step+1 << "/" << this->_max_r_steps
        << ", err in QoI: " << error.l2_norm() << std::endl;
  }

  return;
}

void GRINS::ModelAdaptiveSimulation::print_sim_info()
{
  // Print mesh info if the user wants it
  if( this->_print_mesh_info ) this->_forward_mesh->print_info();

  // Print info if requested
  if( this->_print_equation_system_info )
  {
    out << "\nForward system:\n";
    this->_forward_equation_system->print_info();

    out << "\nAdjoint system:\n";
    this->_adjoint_equation_system->print_info();

    out << "\nResidual system:\n";
    this-> _residual_equation_system->print_info();
  }

  return;
}

std::tr1::shared_ptr<libMesh::EquationSystems> GRINS::ModelAdaptiveSimulation::get_equation_system()
{
  out << "Not implemented for model adaptive sim." << std::endl;
  libmesh_error();
  return _forward_equation_system;
}

void GRINS::ModelAdaptiveSimulation::check_for_restart( const GetPot& input )
{
  const std::string restart_file = input( "restart-options/restart_file", "none" );

  // Most of this was pulled from FIN-S
  if (restart_file != "none")
    {
      std::cout << " ====== Restarting from " << restart_file << std::endl;      

       // Must have correct file type to restart
      if (restart_file.rfind(".xdr") < restart_file.size())
        _forward_equation_system->read(restart_file,libMeshEnums::DECODE,
			      EquationSystems::READ_DATA |
			      EquationSystems::READ_ADDITIONAL_DATA);
      
      else if  (restart_file.rfind(".xda") < restart_file.size())
        _forward_equation_system->read(restart_file,libMeshEnums::READ,
			      EquationSystems::READ_DATA |
			      EquationSystems::READ_ADDITIONAL_DATA);

      else
        {
          std::cerr << "Error: Restart filename must have .xdr or .xda extension!" << std::endl;
          libmesh_error();
        }
      
      GRINS::MultiphysicsSystem& system = 
        _forward_equation_system->get_system<GRINS::MultiphysicsSystem>( "GRINS" );

      // Update the old data
      system.update();
    }

  return;
}

void GRINS::ModelAdaptiveSimulation::attach_neumann_bc_funcs( std::map< std::string, GRINS::NBCContainer > neumann_bcs,
						 GRINS::MultiphysicsSystem* system )
{
  //_neumann_bc_funcs = neumann_bcs;

  if( neumann_bcs.size() > 0 )
    {
      for( std::map< std::string, GRINS::NBCContainer >::iterator bc = neumann_bcs.begin();
	   bc != neumann_bcs.end();
	   bc++ )
	{
	  std::tr1::shared_ptr<GRINS::Physics> physics = system->get_physics( bc->first );
	  physics->attach_neumann_bound_func( bc->second );
	}
    }

  return;
}

void GRINS::ModelAdaptiveSimulation::attach_dirichlet_bc_funcs( std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > dbc_map,
						   GRINS::MultiphysicsSystem* system )
{
  for( std::multimap< GRINS::PhysicsName, GRINS::DBCContainer >::const_iterator it = dbc_map.begin();
       it != dbc_map.end();
       it++ )
    {
      std::tr1::shared_ptr<GRINS::Physics> physics = system->get_physics( it->first );
      
      physics->attach_dirichlet_bound_func( it->second );
    }
  return;
}

#ifdef USE_GRVY_TIMERS
void GRINS::ModelAdaptiveSimulation::attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer )
{
  this->_multiphysics_system->attach_grvy_timer( grvy_timer );
  return;
}
#endif
