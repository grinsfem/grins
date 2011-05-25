//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - a low Mach number Navier-Stokes Finite-Element Solver
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
// Definitions for the MultphysicsSystem class.
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "multiphysics_sys.h"

#include "inc_navier_stokes.h"

GRINS::MultiphysicsSystem::~MultiphysicsSystem()
{
  // Physics* objects get new'ed in the read_input_options call, so we
  // need to delete them.
  for( GRINS::PhysicsListIter physics_iter = _physics_list.begin();
       physics_iter != _physics_list.end();
       physics_iter++ )
    {
      delete (physics_iter->second);
    }
  
  return;
}

void GRINS::MultiphysicsSystem::read_input_options( GetPot& input )
{
  // Read options for MultiphysicsSystem first
  this->verify_analytic_jacobians  = input("linear-nonlinear-solver/verify_analytic_jacobians", 0.0 );

  // Figure out how many physics we are enabling
  int num_physics = input.vector_variable_size("Physics/enabled_physics");

  if( num_physics < 1 )
    {
      std::cerr << "Error: Must enable at least one physics model" << std::endl;
      libmesh_error();
    }

  // Go through and create a physics object for each physics we're enabling
  for( int i = 0; i < num_physics; i++ )
    {
      std::string physics_to_add = input("Physics/enabled_physics", "NULL", i );

      //TODO: Do we want to create an enum list instead for all available physics?
      if( physics_to_add == "IncompressibleNavierStokes" )
	{
	  this->_physics_list.push_back( std::make_pair("IncompressibleNavierStokes", new GRINS::IncompressibleNavierStokes) );
	}
      else
	{
	  std::cerr << "Error: Invalid physics name " << physics_to_add << std::endl;
	  libmesh_error();
	}
    }
  
  // Read the input options for each of the physics
  for( GRINS::PhysicsListIter physics_iter = _physics_list.begin();
       physics_iter != _physics_list.end();
       physics_iter++ )
    {
      (physics_iter->second)->read_input_options( input );
    }

  return;
}

void GRINS::MultiphysicsSystem::init_data()
{
  // First, initalize all the variables. We pass this pointer for the system.
  for( GRINS::PhysicsListIter physics_iter = _physics_list.begin();
       physics_iter != _physics_list.end();
       physics_iter++ )
    {
      (physics_iter->second)->init_variables( this );
    }

  // Next, call register_variable_indices in each physics. We do this
  // by buidling up a "global" variable map for all the physics, then
  // pass this to each physics so it can grab the varible it wants.
  for( GRINS::PhysicsListIter physics_iter = _physics_list.begin();
       physics_iter != _physics_list.end();
       physics_iter++ )
    {
      VariableMap local_map = (physics_iter->second)->get_variable_indices_map();

      _global_map.insert( local_map.begin(), local_map.end() );
    }

  //this->dump_global_variable_map();

  for( GRINS::PhysicsListIter physics_iter = _physics_list.begin();
       physics_iter != _physics_list.end();
       physics_iter++ )
    {
      (physics_iter->second)->register_variable_indices( _global_map );
    }

  // Next, call parent init_data function to intialize everything.
  libMesh::FEMSystem::init_data();

  // Now set time_evolving variables
  for( GRINS::PhysicsListIter physics_iter = _physics_list.begin();
       physics_iter != _physics_list.end();
       physics_iter++ )
    {
      (physics_iter->second)->set_time_evolving_vars( this );
    }

  return;
}

void GRINS::MultiphysicsSystem::init_context( libMesh::DiffContext &context )
{
  //Loop over each physics to initialize relevant variable structures for assembling system
  for( GRINS::PhysicsListIter physics_iter = _physics_list.begin();
       physics_iter != _physics_list.end();
       physics_iter++ )
    {
      (physics_iter->second)->init_context( context );
    }

  return;
}

bool GRINS::MultiphysicsSystem::element_time_derivative( bool request_jacobian,
							 libMesh::DiffContext& context )
{
  // Loop over each physics and compute their contributions
  for( GRINS::PhysicsListIter physics_iter = _physics_list.begin();
       physics_iter != _physics_list.end();
       physics_iter++ )
    {
      (physics_iter->second)->element_time_derivative( request_jacobian, context, this );
    }

  // TODO: Need to think about the implications of this because there might be some
  // TODO: jacobian terms we don't want to compute for efficiency reasons
  return request_jacobian;
}

bool GRINS::MultiphysicsSystem::side_time_derivative( bool request_jacobian,
						      libMesh::DiffContext& context )
{
  // Loop over each physics and compute their contributions
  for( GRINS::PhysicsListIter physics_iter = _physics_list.begin();
       physics_iter != _physics_list.end();
       physics_iter++ )
    {
      (physics_iter->second)->side_time_derivative( request_jacobian, context, this );
    }

  // TODO: Need to think about the implications of this because there might be some
  // TODO: jacobian terms we don't want to compute for efficiency reasons
  return request_jacobian;
}

bool GRINS::MultiphysicsSystem::element_constraint( bool request_jacobian,
						    libMesh::DiffContext& context )
{
  // Loop over each physics and compute their contributions
  for( GRINS::PhysicsListIter physics_iter = _physics_list.begin();
       physics_iter != _physics_list.end();
       physics_iter++ )
    {
      (physics_iter->second)->element_constraint( request_jacobian, context, this );
    }

  // TODO: Need to think about the implications of this because there might be some
  // TODO: jacobian terms we don't want to compute for efficiency reasons
  return request_jacobian;
}

bool GRINS::MultiphysicsSystem::side_constraint( bool request_jacobian,
						 libMesh::DiffContext& context )
{
  // Loop over each physics and compute their contributions
  for( GRINS::PhysicsListIter physics_iter = _physics_list.begin();
       physics_iter != _physics_list.end();
       physics_iter++ )
    {
      (physics_iter->second)->side_constraint( request_jacobian, context, this );
    }

  // TODO: Need to think about the implications of this because there might be some
  // TODO: jacobian terms we don't want to compute for efficiency reasons
  return request_jacobian;
}

bool GRINS::MultiphysicsSystem::mass_residual( bool request_jacobian,
					       libMesh::DiffContext& context )
{
  // Loop over each physics and compute their contributions
  for( GRINS::PhysicsListIter physics_iter = _physics_list.begin();
       physics_iter != _physics_list.end();
       physics_iter++ )
    {
      (physics_iter->second)->mass_residual( request_jacobian, context, this );
    }

  // TODO: Need to think about the implications of this because there might be some
  // TODO: jacobian terms we don't want to compute for efficiency reasons
  return request_jacobian;
}

void GRINS::MultiphysicsSystem::dump_global_variable_map( )
{
  std::cout << "Dumping Variable map" << std::endl;
  for( GRINS::VariableMapConstIt it = _global_map.begin();
       it != _global_map.end();
       it++ )
    {
      std::cout << "Variable: " << it->first
		<< ", Index: " << it->second << std::endl; 
    }
  std::cout << "Done dumping variable map" << std::endl;

  return;
}

#ifdef USE_GRVY_TIMERS
void GRINS::MultiphysicsSystem::attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer )
{
  _timer = grvy_timer;

  // Attach timers to each physics
  for( GRINS::PhysicsListIter physics_iter = _physics_list.begin();
       physics_iter != _physics_list.end();
       physics_iter++ )
    {
      (physics_iter->second)->attach_grvy_timer( grvy_timer );
    }

  return;
}
#endif
