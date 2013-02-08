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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "grins/physics.h"

// GRINS
#include "grins/bc_handling_base.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/elem.h"

namespace GRINS
{
  // Initialize static members
  bool Physics::_is_steady = false;

  Physics::Physics( const std::string& physics_name,
		    const GetPot& input )
    : _physics_name( physics_name ),
      _bc_handler(NULL)
  {
    this->read_input_options(input);
    return;
  }

  Physics::~Physics()
  {
    // If a derived class created a bc_handler object, we kill it here.
    if( _bc_handler ) delete _bc_handler;
    return;
  }

  void Physics::read_input_options( const GetPot& input )
  {
    int num_ids = input.vector_variable_size( "Physics/"+this->_physics_name+"/enabled_subdomains" );

    for( int i = 0; i < num_ids; i++ )
      {
	libMesh::subdomain_id_type dumvar = input( "Physics/"+this->_physics_name+"/enabled_subdomains", -1, i );
	_enabled_subdomains.insert( dumvar );
      }

    return;
  }

  bool Physics::enabled_on_elem( const libMesh::Elem* elem )
  {
    // Check if enabled_subdomains flag has been set
    if( _enabled_subdomains.empty() ) 
      return true;
  
    // Check if current physics is enabled on elem
    if( _enabled_subdomains.find( elem->subdomain_id() ) == _enabled_subdomains.end() )
      return false;

    return true;
  }

  void Physics::set_is_steady( bool is_steady )
  {
    _is_steady = is_steady;
    return;
  }

  bool Physics::is_steady() const
  {
    return _is_steady;
  }

  void Physics::set_time_evolving_vars( libMesh::FEMSystem* /*system*/ )
  {
    return;
  }

  void Physics::init_bcs( libMesh::FEMSystem* system )
  {
    // Only need to init BC's if the physics actually created a handler
    if( _bc_handler )
      {
	_bc_handler->init_bc_data( *system );
	_bc_handler->init_dirichlet_bcs( system );
	_bc_handler->init_dirichlet_bc_func_objs( system );
	_bc_handler->init_periodic_bcs( system );
      }

    return;
  }

  void Physics::attach_neumann_bound_func( NBCContainer& neumann_bcs )
  {
    _bc_handler->attach_neumann_bound_func( neumann_bcs );
    return;
  }

  void Physics::attach_dirichlet_bound_func( const DBCContainer& dirichlet_bc )
  {
    _bc_handler->attach_dirichlet_bound_func( dirichlet_bc );
    return;
  }
  
  void Physics::init_context( libMesh::FEMContext& /*context*/ )
  {
    return;
  }

  void Physics::compute_element_time_derivative_cache( const libMesh::FEMContext&,
						       CachedValues& ) const
  {
    return;
  }

  void Physics::compute_side_time_derivative_cache( const libMesh::FEMContext& /*context*/,
						    CachedValues& /*cache*/ ) const
   {
     return;
   }

  void Physics::compute_element_constraint_cache( const libMesh::FEMContext& /*context*/,
						  CachedValues& /*cache*/ ) const
  {
    return;
  }

  void Physics::compute_side_constraint_cache( const libMesh::FEMContext& /*context*/,
					       CachedValues& /*cache*/ ) const
  {
    return;
  }

  void Physics::compute_mass_residual_cache( const libMesh::FEMContext& /*context*/,
					     CachedValues& /*cache*/ ) const
  {
    return;
  }

  void Physics::compute_element_cache( const libMesh::FEMContext&,
				       const std::vector<libMesh::Point>&,
				       CachedValues& ) const
  {
    return;
  }

  void Physics::element_time_derivative( bool /*compute_jacobian*/,
					 libMesh::FEMContext& /*context*/,
					 CachedValues& /*cache*/ )
  {
    return;
  }

  void Physics::side_time_derivative( bool /*compute_jacobian*/,
				      libMesh::FEMContext& /*context*/,
				      CachedValues& /*cache*/ )
  {
    return;
  }

  void Physics::element_constraint( bool /*compute_jacobian*/,
				    libMesh::FEMContext& /*context*/,
				    CachedValues& /*cache*/ )
  {
    return;
  }

  void Physics::side_constraint( bool /*compute_jacobian*/,
				 libMesh::FEMContext& /*context*/,
				 CachedValues& /*cache*/ )
  {
    return;
  }   

  void Physics::mass_residual( bool /*compute_jacobian*/,
			       libMesh::FEMContext& /*context*/,
			       CachedValues& /*cache*/ )
  {
    return;
  }

#ifdef GRINS_USE_GRVY_TIMERS
  void Physics::attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer )
  {
    _timer = grvy_timer;
    return;
  }
#endif

} // namespace GRINS
