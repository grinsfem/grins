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
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{

  MultiphysicsSystem::MultiphysicsSystem( libMesh::EquationSystems& es,
					  const std::string& name,
					  const unsigned int number )
    : FEMSystem(es, name, number),
      _use_numerical_jacobians_only(false)
  {
    return;
  }

  MultiphysicsSystem::~MultiphysicsSystem()
  {
    return;
  }

  void MultiphysicsSystem::attach_physics_list( PhysicsList physics_list )
  {
    _physics_list = physics_list;
    return;
  }

  void MultiphysicsSystem::read_input_options( const GetPot& input )
  {
    // Read options for MultiphysicsSystem first
    this->verify_analytic_jacobians = input("linear-nonlinear-solver/verify_analytic_jacobians", 0.0 );
    this->print_element_jacobians = input("screen-options/print_element_jacobians", false );
    _use_numerical_jacobians_only = input("linear-nonlinear-solver/use_numerical_jacobians_only", false );
  }

  void MultiphysicsSystem::init_data()
  {
    // Need this to be true because of our overloading of the
    // mass_residual function.
    // This is data in FEMSystem. MUST be set before FEMSystem::init_data.
    use_fixed_solution = true;

    // Initalize all the variables. We pass this pointer for the system.
    /* NOTE: We CANNOT fuse this loop with the others. This loop
       MUST complete first. */
    /*! \todo Figure out how to tell compilers not to fuse this loop when
      they want to be aggressive. */
    for( PhysicsListIter physics_iter = _physics_list.begin();
	 physics_iter != _physics_list.end();
	 physics_iter++ )
      {
	(physics_iter->second)->init_variables( this );
      }

    // Now set time_evolving variables
    for( PhysicsListIter physics_iter = _physics_list.begin();
	 physics_iter != _physics_list.end();
	 physics_iter++ )
      {
	(physics_iter->second)->set_time_evolving_vars( this );
      }

    // Set whether the problem we're solving is steady or not
    // Since the variable is static, just call one Physics class
    {
      (_physics_list.begin()->second)->set_is_steady((this->time_solver)->is_steady());
    }

    for( PhysicsListIter physics_iter = _physics_list.begin();
	 physics_iter != _physics_list.end();
	 physics_iter++ )
      {
	// Initialize builtin BC's for each physics
	(physics_iter->second)->init_bcs( this );
      }

    // Next, call parent init_data function to intialize everything.
    libMesh::FEMSystem::init_data();

    return;
  }

  void MultiphysicsSystem::init_context( libMesh::DiffContext &context )
  {
    libMesh::FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

    //Loop over each physics to initialize relevant variable structures for assembling system
    for( PhysicsListIter physics_iter = _physics_list.begin();
	 physics_iter != _physics_list.end();
	 physics_iter++ )
      {
	(physics_iter->second)->init_context( c );
      }

    return;
  }

  bool MultiphysicsSystem::element_time_derivative( bool request_jacobian,
						    libMesh::DiffContext& context )
  {
    libMesh::FEMContext& c = libmesh_cast_ref<libMesh::FEMContext&>( context );
  
    bool compute_jacobian = true;
    if( !request_jacobian || _use_numerical_jacobians_only ) compute_jacobian = false;

    CachedValues cache;

    // Now compute cache for this element
    for( PhysicsListIter physics_iter = _physics_list.begin();
	 physics_iter != _physics_list.end();
	 physics_iter++ )
      {
	(physics_iter->second)->compute_element_time_derivative_cache( c, cache );
      }

    // Loop over each physics and compute their contributions
    for( PhysicsListIter physics_iter = _physics_list.begin();
	 physics_iter != _physics_list.end();
	 physics_iter++ )
      {
	// Only compute if physics is active on current subdomain or globally
	if( (physics_iter->second)->enabled_on_elem( c.elem ) )
	  {
	    (physics_iter->second)->element_time_derivative( compute_jacobian, c,
							     cache );
	  }
      }

    // TODO: Need to think about the implications of this because there might be some
    // TODO: jacobian terms we don't want to compute for efficiency reasons
    return compute_jacobian;
  }

  bool MultiphysicsSystem::side_time_derivative( bool request_jacobian,
						 libMesh::DiffContext& context )
  {
    libMesh::FEMContext& c = libmesh_cast_ref<libMesh::FEMContext&>( context );

    bool compute_jacobian = true;
    if( !request_jacobian || _use_numerical_jacobians_only ) compute_jacobian = false;

    CachedValues cache;

    // Now compute cache for this element
    for( PhysicsListIter physics_iter = _physics_list.begin();
	 physics_iter != _physics_list.end();
	 physics_iter++ )
      {
	(physics_iter->second)->compute_side_time_derivative_cache( c, cache );
      }

    // Loop over each physics and compute their contributions
    for( PhysicsListIter physics_iter = _physics_list.begin();
	 physics_iter != _physics_list.end();
	 physics_iter++ )
      {
	// Only compute if physics is active on current subdomain or globally
	if( (physics_iter->second)->enabled_on_elem( c.elem ) )
	  {
	    (physics_iter->second)->side_time_derivative( compute_jacobian, c,
							  cache );
	  }
      }

    // TODO: Need to think about the implications of this because there might be some
    // TODO: jacobian terms we don't want to compute for efficiency reasons
    return compute_jacobian;
  }

  bool MultiphysicsSystem::element_constraint( bool request_jacobian,
					       libMesh::DiffContext& context )
  {
    libMesh::FEMContext& c = libmesh_cast_ref<libMesh::FEMContext&>( context );

    bool compute_jacobian = true;
    if( !request_jacobian || _use_numerical_jacobians_only ) compute_jacobian = false;

    CachedValues cache;

    // Now compute cache for this element
    for( PhysicsListIter physics_iter = _physics_list.begin();
	 physics_iter != _physics_list.end();
	 physics_iter++ )
      {
	(physics_iter->second)->compute_element_constraint_cache( c, cache );
      }

    // Loop over each physics and compute their contributions
    for( PhysicsListIter physics_iter = _physics_list.begin();
	 physics_iter != _physics_list.end();
	 physics_iter++ )
      {
	// Only compute if physics is active on current subdomain or globally
	if( (physics_iter->second)->enabled_on_elem( c.elem ) )
	  {
	    (physics_iter->second)->element_constraint( compute_jacobian, c,
							cache);
	  }
      }

    // TODO: Need to think about the implications of this because there might be some
    // TODO: jacobian terms we don't want to compute for efficiency reasons
    return compute_jacobian;
  }

  bool MultiphysicsSystem::side_constraint( bool request_jacobian,
					    libMesh::DiffContext& context )
  {
    libMesh::FEMContext& c = libmesh_cast_ref<libMesh::FEMContext&>( context );

    bool compute_jacobian = true;
    if( !request_jacobian || _use_numerical_jacobians_only ) compute_jacobian = false;

    CachedValues cache;

    // Now compute cache for this element
    for( PhysicsListIter physics_iter = _physics_list.begin();
	 physics_iter != _physics_list.end();
	 physics_iter++ )
      {
	(physics_iter->second)->compute_side_constraint_cache( c, cache );
      }

    // Loop over each physics and compute their contributions
    for( PhysicsListIter physics_iter = _physics_list.begin();
	 physics_iter != _physics_list.end();
	 physics_iter++ )
      {
	// Only compute if physics is active on current subdomain or globally
	if( (physics_iter->second)->enabled_on_elem( c.elem ) )
	  {
	    (physics_iter->second)->side_constraint( compute_jacobian, c,
						     cache);
	  }
      }

    // TODO: Need to think about the implications of this because there might be some
    // TODO: jacobian terms we don't want to compute for efficiency reasons
    return compute_jacobian;
  }

  bool MultiphysicsSystem::mass_residual( bool request_jacobian,
					  libMesh::DiffContext& context )
  {
    libMesh::FEMContext& c = libmesh_cast_ref<libMesh::FEMContext&>( context );

    bool compute_jacobian = true;
    if( !request_jacobian || _use_numerical_jacobians_only ) compute_jacobian = false;

    CachedValues cache;

    // Now compute cache for this element
    for( PhysicsListIter physics_iter = _physics_list.begin();
	 physics_iter != _physics_list.end();
	 physics_iter++ )
      {
	(physics_iter->second)->compute_mass_residual_cache( c, cache );
      }

    // Loop over each physics and compute their contributions
    for( PhysicsListIter physics_iter = _physics_list.begin();
	 physics_iter != _physics_list.end();
	 physics_iter++ )
      {
	// Only compute if physics is active on current subdomain or globally
	if( (physics_iter->second)->enabled_on_elem( c.elem ) )
	  {
	    (physics_iter->second)->mass_residual( compute_jacobian, c,
						   cache);
	  }
      }

    // TODO: Need to think about the implications of this because there might be some
    // TODO: jacobian terms we don't want to compute for efficiency reasons
    return compute_jacobian;
  }

  std::tr1::shared_ptr<Physics> MultiphysicsSystem::get_physics( const std::string physics_name )
  {
    if( _physics_list.find( physics_name ) == _physics_list.end() )
      {
	std::cerr << "Error: Could not find physics " << physics_name << std::endl;
	libmesh_error();
      }

    return _physics_list[physics_name];
  }

  bool MultiphysicsSystem::has_physics( const std::string physics_name ) const
  {
    bool has_physics = false;

    if( _physics_list.find(physics_name) != _physics_list.end() )
      has_physics = true;

    return has_physics;
  }

  void MultiphysicsSystem::compute_element_cache( const libMesh::FEMContext& context,
						  const std::vector<libMesh::Point>& points,
						  CachedValues& cache ) const
  {
    for( PhysicsListIter physics_iter = _physics_list.begin();
	 physics_iter != _physics_list.end();
	 physics_iter++ )
      {
	(physics_iter->second)->compute_element_cache( context, points, cache );
      }
    return;
  }

#ifdef GRINS_USE_GRVY_TIMERS
  void MultiphysicsSystem::attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer )
  {
    _timer = grvy_timer;

    // Attach timers to each physics
    for( PhysicsListIter physics_iter = _physics_list.begin();
	 physics_iter != _physics_list.end();
	 physics_iter++ )
      {
	(physics_iter->second)->attach_grvy_timer( grvy_timer );
      }

    return;
  }
#endif


} // namespace GRINS
