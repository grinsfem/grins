//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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
#include "grins/ic_handling_base.h"
#include "grins/fe_variables_base.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"

namespace GRINS
{
  // Initialize static members
  bool Physics::_is_steady = false;
  bool Physics::_is_axisymmetric = false;

  Physics::Physics( const std::string& physics_name,
                    const GetPot& input )
    : ParameterUser(physics_name),
      _physics_name( physics_name ),
      _ic_handler(new ICHandlingBase(physics_name))
  {
    this->parse_enabled_subdomains(input,physics_name);

    // Check if this is an axisymmetric problem
    // There will be redundant calls for multiple Physics objects,
    // but this guarantees we've parsed is_axisymmetric and then
    // we can check in subclasses and error out if the Physics
    // doesn't support axisymmetry.
    Physics::set_is_axisymmetric( input("Physics/is_axisymmetric",false) );
  }

  Physics::~Physics()
  {
    if( _ic_handler ) delete _ic_handler;
  }

  void Physics::parse_enabled_subdomains( const GetPot& input,
                                          const std::string& physics_name )
  {
    int num_ids = input.vector_variable_size( "Physics/"+physics_name+"/enabled_subdomains" );

    for( int i = 0; i < num_ids; i++ )
      {
        libMesh::subdomain_id_type dumvar = input( "Physics/"+physics_name+"/enabled_subdomains", -1, i );
        _enabled_subdomains.insert( dumvar );
      }
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


  void Physics::init_ics( libMesh::FEMSystem* system,
                          libMesh::CompositeFunction<libMesh::Number>& all_ics )
  {
    if( _ic_handler )
      {
        _ic_handler->init_ic_data( *system, all_ics );
      }

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

  void Physics::compute_postprocessed_quantity( unsigned int /*quantity_index*/,
                                                const AssemblyContext& /*context*/,
                                                const libMesh::Point& /*point*/,
                                                libMesh::Real& /*value*/ )
  {
    return;
  }

  std::unique_ptr<libMesh::FEGenericBase<libMesh::Real> > Physics::build_new_fe( const libMesh::Elem* elem,
                                                                                 const libMesh::FEGenericBase<libMesh::Real>* fe,
                                                                                 const libMesh::Point p )
  {
    using namespace libMesh;
    FEType fe_type = fe->get_fe_type();

    // If we don't have an Elem to evaluate on, then the only functions
    // we can sensibly evaluate are the scalar dofs which are the same
    // everywhere.
    libmesh_assert(elem || fe_type.family == SCALAR);

    unsigned int elem_dim = elem ? elem->dim() : 0;

    std::unique_ptr<FEGenericBase<libMesh::Real> >
      fe_new(FEGenericBase<libMesh::Real>::build(elem_dim, fe_type));

    // Map the physical co-ordinates to the master co-ordinates using the inverse_map from fe_interface.h
    // Build a vector of point co-ordinates to send to reinit
    Point master_point = elem ?
      FEInterface::inverse_map(elem_dim, fe_type, elem, p) :
      Point(0);

    std::vector<Point> coor(1, master_point);

    // Reinitialize the element and compute the shape function values at coor
    fe_new->reinit (elem, &coor);

    return fe_new;
  }

  void Physics::check_var_subdomain_consistency( const FEVariablesBase& var ) const
  {
    const std::set<libMesh::subdomain_id_type>& var_subdomains = var.subdomain_ids();

    // If both are empty or only the var is empty, we don't need to do anything
    // since the check is automatically satisified. Empty in both cases means
    // enabled on all subdomains

    if( _enabled_subdomains.empty() && !var_subdomains.empty() )
      libmesh_error_msg("ERROR: Physics enabled on all subdomains but variable is not!");

    if( !_enabled_subdomains.empty() && !var_subdomains.empty() )
      for( std::set<libMesh::subdomain_id_type>::const_iterator it = _enabled_subdomains.begin(); it != _enabled_subdomains.end(); ++it )
        if( var_subdomains.find(*it) == var_subdomains.end() )
          libmesh_error_msg("ERROR: Could not find subdomain " << *it << " in varaible!");
  }

} // namespace GRINS
