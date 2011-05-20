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
	  this->_physics_list.push_back( new GRINS::IncompressibleNavierStokes );
	}
      else
	{
	  std::cerr << "Error: Invalid physics name " << physics_to_add << std::endl;
	  libmesh_error();
	}
    }
  
  // Read the input options for each of the physics
  for( std::vector<Physics*>::iterator physics = _physics_list.begin();
       physics != _physics_list.end();
       physics++ )
    {
      (*physics)->read_input_options( input );
    }

  return;
}

void GRINS::MultiphysicsSystem::init_data()
{

  // First, initalize all the variables. We pass this pointer for the system.
  for( std::vector<Physics*>::iterator physics = _physics_list.begin();
       physics != _physics_list.end();
       physics++ )
    {
      (*physics)->init_variables( this );
    }

  // Next, call parent init_data function to intialize everything.
  libMesh::FEMSystem::init_data();

  // Now set time_evolving variables
  for( std::vector<Physics*>::iterator physics = _physics_list.begin();
       physics != _physics_list.end();
       physics++ )
    {
      (*physics)->set_time_evolving_vars( this );
    }

  return;
}

void GRINS::MultiphysicsSystem::init_context( libMesh::DiffContext &context )
{
  //Loop over each physics to initialize relevant variable structures for assembling system
  for( std::vector<Physics*>::iterator physics = _physics_list.begin();
       physics != _physics_list.end();
       physics++ )
    {
      (*physics)->init_context( context );
    }

  return;
}

bool GRINS::MultiphysicsSystem::element_time_derivative( bool request_jacobian,
							 libMesh::DiffContext& context )
{
  // Loop over each physics and compute their contributions
  for( std::vector<Physics*>::iterator physics = _physics_list.begin();
       physics != _physics_list.end();
       physics++ )
    {
      (*physics)->element_time_derivative( request_jacobian, context, this );
    }

  // TODO: Need to think about the implications of this because there might be some
  // TODO: jacobian terms we don't want to compute for efficiency reasons
  return request_jacobian;
}

bool GRINS::MultiphysicsSystem::side_time_derivative( bool request_jacobian,
						      libMesh::DiffContext& context )
{
  // Loop over each physics and compute their contributions
  for( std::vector<Physics*>::iterator physics = _physics_list.begin();
       physics != _physics_list.end();
       physics++ )
    {
      (*physics)->side_time_derivative( request_jacobian, context, this );
    }

  // TODO: Need to think about the implications of this because there might be some
  // TODO: jacobian terms we don't want to compute for efficiency reasons
  return request_jacobian;
}

bool GRINS::MultiphysicsSystem::element_constraint( bool request_jacobian,
						    libMesh::DiffContext& context )
{
  // Loop over each physics and compute their contributions
  for( std::vector<Physics*>::iterator physics = _physics_list.begin();
       physics != _physics_list.end();
       physics++ )
    {
      (*physics)->element_constraint( request_jacobian, context, this );
    }

  // TODO: Need to think about the implications of this because there might be some
  // TODO: jacobian terms we don't want to compute for efficiency reasons
  return request_jacobian;
}

bool GRINS::MultiphysicsSystem::side_constraint( bool request_jacobian,
						 libMesh::DiffContext& context )
{
  // Loop over each physics and compute their contributions
  for( std::vector<Physics*>::iterator physics = _physics_list.begin();
       physics != _physics_list.end();
       physics++ )
    {
      (*physics)->side_constraint( request_jacobian, context, this );
    }

  // TODO: Need to think about the implications of this because there might be some
  // TODO: jacobian terms we don't want to compute for efficiency reasons
  return request_jacobian;
}

bool GRINS::MultiphysicsSystem::mass_residual( bool request_jacobian,
					       libMesh::DiffContext& context )
{
  // Loop over each physics and compute their contributions
  for( std::vector<Physics*>::iterator physics = _physics_list.begin();
       physics != _physics_list.end();
       physics++ )
    {
      (*physics)->mass_residual( request_jacobian, context, this );
    }

  // TODO: Need to think about the implications of this because there might be some
  // TODO: jacobian terms we don't want to compute for efficiency reasons
  return request_jacobian;
}

