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

#include "multiphysics_sys.h"

GRINS::MultiphysicsSystem::MultiphysicsSystem( libMesh::EquationSystems& es,
					       const std::string& name,
					       const unsigned int number )
  : FEMSystem(es, name, number)
    {}

GRINS::MultiphysicsSystem::~MultiphysicsSystem()
{
  return;
}

void GRINS::MultiphysicsSystem::attach_physics_list( PhysicsList physics_list )
{
  _physics_list = physics_list;
  return;
}

void GRINS::MultiphysicsSystem::read_input_options( const GetPot& input )
{
  // Read options for MultiphysicsSystem first
  this->verify_analytic_jacobians  = input("linear-nonlinear-solver/verify_analytic_jacobians", 0.0 );

  // Read the input options for each of the physics
  for( GRINS::PhysicsListIter physics_iter = _physics_list.begin();
       physics_iter != _physics_list.end();
       physics_iter++ )
    {
      (physics_iter->second)->read_input_options( input );
    }

  // Read boundary condition data
  for( GRINS::PhysicsListIter physics_iter = _physics_list.begin();
       physics_iter != _physics_list.end();
       physics_iter++ )
    {
      (physics_iter->second)->read_bc_data( input );
    }
  return;
}

void GRINS::MultiphysicsSystem::init_data()
{
  // Need this to be true because of our overloading of the
  // mass_residual function.
  // This is data in FEMSystem. MUST be set before FEMSystem::init_data.
  use_fixed_solution = true;

  // Initalize all the variables. We pass this pointer for the system.
  for( GRINS::PhysicsListIter physics_iter = _physics_list.begin();
       physics_iter != _physics_list.end();
       physics_iter++ )
    {
      (physics_iter->second)->init_variables( this );
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

std::tr1::shared_ptr<GRINS::Physics> GRINS::MultiphysicsSystem::get_physics( const std::string physics_name )
{
  return _physics_list[physics_name];
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
