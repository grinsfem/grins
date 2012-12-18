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

#include "physics.h"

namespace GRINS
{

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


  void Physics::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    return;
  }

  void Physics::init_bcs( libMesh::FEMSystem* system )
  {
    // Only need to init BC's if the physics actually created a handler
    if( _bc_handler )
      {
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

  void Physics::init_cache( CachedValues& )
  {
    return;
  }

  void Physics::compute_cache( libMesh::FEMContext&, CachedValues& )
  {
    return;
  }

  void Physics::compute_cache( libMesh::FEMContext&, CachedValues&,
			       const std::vector<libMesh::Point>& )
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
