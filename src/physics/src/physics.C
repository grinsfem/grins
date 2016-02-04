//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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
#include "grins/physics.h"

// GRINS
#include "grins/bc_handling_base.h"
#include "grins/ic_handling_base.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/elem.h"

namespace GRINS
{
  // Initialize static members
  bool Physics::_is_steady = false;

  Physics::Physics( const std::string& physics_name,
		    const GetPot& input )
    : ParameterUser(physics_name),
      _physics_name( physics_name ),
      _bc_handler(NULL),
      _ic_handler(new ICHandlingBase(physics_name)),
      _is_axisymmetric(false)
  {
    this->read_input_options(input);

    if( input( "Physics/is_axisymmetric", false ) )
      {
        _is_axisymmetric = true;
      }

    return;
  }

  Physics::~Physics()
  {
    // If a derived class created a bc_handler object, we kill it here.
    if( _bc_handler ) delete _bc_handler;

    if( _ic_handler ) delete _ic_handler;

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
    // Check if enabled_subdomains flag has been set and if we're
    // looking at a real element (rather than a nonlocal evaluation)
    if( !elem || _enabled_subdomains.empty() ) 
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

  void Physics::auxiliary_init( MultiphysicsSystem& /*system*/ )
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


  void Physics::init_ics( libMesh::FEMSystem* system,
                          libMesh::CompositeFunction<libMesh::Number>& all_ics )
  {
    if( _ic_handler )
      {
	_ic_handler->init_ic_data( *system, all_ics );
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
  
  void Physics::init_context( AssemblyContext& /*context*/ )
  {
    return;
  }

  void Physics::register_postprocessing_vars( const GetPot& /*input*/,
                                              PostProcessedQuantities<libMesh::Real>& /*postprocessing*/ )
  {
    return;
  }

  void Physics::compute_element_time_derivative_cache( const AssemblyContext&,
						       CachedValues& )
  {
    return;
  }

  void Physics::compute_side_time_derivative_cache( const AssemblyContext& /*context*/,
						    CachedValues& /*cache*/ )
  {
    return;
  }

  void Physics::compute_nonlocal_time_derivative_cache( const AssemblyContext& /*context*/,
						        CachedValues& /*cache*/ )
  {
    return;
  }

  void Physics::compute_element_constraint_cache( const AssemblyContext& /*context*/,
						  CachedValues& /*cache*/ )
  {
    return;
  }

  void Physics::compute_side_constraint_cache( const AssemblyContext& /*context*/,
					       CachedValues& /*cache*/ )
  {
    return;
  }

  void Physics::compute_nonlocal_constraint_cache( const AssemblyContext& /*context*/,
					           CachedValues& /*cache*/ )
  {
    return;
  }

  void Physics::compute_damping_residual_cache( const AssemblyContext& /*context*/,
                                                CachedValues& /*cache*/ )
  {
    return;
  }

  void Physics::compute_mass_residual_cache( const AssemblyContext& /*context*/,
					     CachedValues& /*cache*/ )
  {
    return;
  }

  void Physics::compute_nonlocal_mass_residual_cache( const AssemblyContext& /*context*/,
					              CachedValues& /*cache*/ )
  {
    return;
  }

  void Physics::element_time_derivative( bool /*compute_jacobian*/,
					 AssemblyContext& /*context*/,
					 CachedValues& /*cache*/ )
  {
    return;
  }

  void Physics::side_time_derivative( bool /*compute_jacobian*/,
				      AssemblyContext& /*context*/,
				      CachedValues& /*cache*/ )
  {
    return;
  }

  void Physics::nonlocal_time_derivative( bool /*compute_jacobian*/,
				          AssemblyContext& /*context*/,
				          CachedValues& /*cache*/ )
  {
    return;
  }

  void Physics::element_constraint( bool /*compute_jacobian*/,
				    AssemblyContext& /*context*/,
				    CachedValues& /*cache*/ )
  {
    return;
  }

  void Physics::side_constraint( bool /*compute_jacobian*/,
				 AssemblyContext& /*context*/,
				 CachedValues& /*cache*/ )
  {
    return;
  }   

  void Physics::nonlocal_constraint( bool /*compute_jacobian*/,
				     AssemblyContext& /*context*/,
				     CachedValues& /*cache*/ )
  {
    return;
  }

  void Physics::damping_residual( bool /*compute_jacobian*/,
                                  AssemblyContext& /*context*/,
                                  CachedValues& /*cache*/ )
  {
    return;
  }

  void Physics::mass_residual( bool /*compute_jacobian*/,
			       AssemblyContext& /*context*/,
			       CachedValues& /*cache*/ )
  {
    return;
  }

  void Physics::nonlocal_mass_residual( bool /*compute_jacobian*/,
			                AssemblyContext& /*context*/,
			                CachedValues& /*cache*/ )
  {
    return;
  }

  void Physics::compute_postprocessed_quantity( unsigned int /*quantity_index*/,
                                                const AssemblyContext& /*context*/,
                                                const libMesh::Point& /*point*/,
                                                libMesh::Real& /*value*/ )
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
