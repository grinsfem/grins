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
#include "grins/bc_handling_base.h"

// libMesh
#include "libmesh/fem_context.h"
#include "libmesh/fem_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/dof_map.h"

namespace GRINS
{
  bool BCHandlingBase::_axisymmetric = false; 

  BCHandlingBase::BCHandlingBase(const std::string& physics_name)
    : _num_periodic_bcs(0),
      _physics_name( physics_name )   
  {
    return;
  }

  BCHandlingBase::~BCHandlingBase()
  {
    return;
  }

  void BCHandlingBase::attach_neumann_bound_func( NBCContainer& neumann_bcs )
  {
    _neumann_bound_funcs.insert( std::make_pair( neumann_bcs.get_bc_id(), neumann_bcs ) );
    return;
  }

  void BCHandlingBase::attach_dirichlet_bound_func( const DBCContainer& dirichlet_bc )
  {
    _dirichlet_bound_funcs.push_back( dirichlet_bc );
    return;
  }

  void BCHandlingBase::read_bc_data( const GetPot& input, const std::string& id_str,
				     const std::string& bc_str)
  {
    int num_ids = input.vector_variable_size(id_str);
    int num_bcs = input.vector_variable_size(bc_str);

    if( num_ids != num_bcs )
      {
	std::cerr << "Error: Must specify equal number of boundary ids and boundary conditions"
		  << std::endl;
	libmesh_error();
      }

    for( int i = 0; i < num_ids; i++ )
      {
	int bc_id = input(id_str, -1, i );
	std::string bc_type_in = input(bc_str, "NULL", i );

	int bc_type = this->string_to_int( bc_type_in );

	std::stringstream ss;
	ss << bc_id;
	std::string bc_id_string = ss.str();

	this->init_bc_types( bc_id, bc_id_string, bc_type, input );
      }

    return;
  }

  void BCHandlingBase::init_bc_data( const libMesh::FEMSystem& /*system*/ )
  {
    return;
  }

  void BCHandlingBase::init_dirichlet_bc_func_objs( libMesh::FEMSystem* system ) const
  {
    libMesh::DofMap& dof_map = system->get_dof_map();

    for( std::vector< DBCContainer >::const_iterator 
	   it = _dirichlet_bound_funcs.begin();
	 it != _dirichlet_bound_funcs.end();
	 it++ )
      {
	// First, get variable names and convert them to variable id's
	std::vector<VariableName> var_names = (*it).get_var_names();
      
	std::vector<VariableIndex> dbc_vars;

	for( std::vector<VariableName>::const_iterator name = var_names.begin();
	     name != var_names.end();
	     name++ )
	  {
	    dbc_vars.push_back( system->variable_number( *name ) );
	  }
      
	// Get bc_id's
	std::set<BoundaryID> bc_ids = (*it).get_bc_ids();
      
	// Get Dirichlet bc functor
	std::tr1::shared_ptr<libMesh::FunctionBase<Number> > func = (*it).get_func();

	// Now create DirichletBoundary object and give it to libMesh
	// libMesh makes it own copy of the DirichletBoundary so we can
	// let this one die.
	libMesh::DirichletBoundary dbc( bc_ids, dbc_vars, &*func );
      
	dof_map.add_dirichlet_boundary( dbc );
      }

    return;
  }

  void BCHandlingBase::init_dirichlet_bcs( libMesh::FEMSystem* system ) const
  {
    libMesh::DofMap& dof_map = system->get_dof_map();

    for( std::map< BoundaryID,BCType >::const_iterator it = _dirichlet_bc_map.begin();
	 it != _dirichlet_bc_map.end();
	 it++ )
      {
	this->user_init_dirichlet_bcs( system, dof_map, it->first, it->second );
      }

    return;
  }

  void BCHandlingBase::init_periodic_bcs( libMesh::FEMSystem* system ) const
  {
    /* Consistency check required to make sure user actually set periodic bc data.
       This is needed because of the way we parse periodic boundary conditions. */
    if( _periodic_bcs.size() != _num_periodic_bcs/2 )
      {
	std::cerr << "=========================================================="  << std::endl
		  << "Error: Inconsistency in perioidic boundary conditions."      << std::endl
		  << "       There were " << _num_periodic_bcs << " deteced but "  << std::endl
		  << "       only " << _periodic_bcs.size() << " pairs of data "   << std::endl
		  << "       were set."                                            << std::endl
		  << "=========================================================="  << std::endl;
	libmesh_error();
      }

    libMesh::DofMap& dof_map = system->get_dof_map();

    for( std::vector< PBCContainer >::const_iterator it = _periodic_bcs.begin();
	 it != _periodic_bcs.end();
	 it++ )
      {
	libMesh::PeriodicBoundary bc( it->get_offset_vector() );
	bc.myboundary = it->get_master_bcid();
	bc.pairedboundary = it->get_slave_bcid();

	dof_map.add_periodic_boundary( bc );
      }

    return;
  }

  void BCHandlingBase::apply_neumann_bcs( libMesh::FEMContext& context,
					  const GRINS::CachedValues& cache,
					  const bool request_jacobian,
					  const BoundaryID bc_id ) const
  {
    std::map< BoundaryID, BCType>::const_iterator 
      bc_map_it = _neumann_bc_map.find( bc_id );

    /* We assume that if you didn't put a boundary id in, then you didn't want to
       set a boundary condition on that boundary. */
    if( bc_map_it != _neumann_bc_map.end() )
      {
	this->user_apply_neumann_bcs( context, cache, request_jacobian,
				      bc_id, bc_map_it->second );
      }
    return;
  }

  void BCHandlingBase::set_dirichlet_bc_type( BoundaryID bc_id, int bc_type )
  {
    _dirichlet_bc_map[bc_id] = bc_type;
    return;
  }

  void BCHandlingBase::set_neumann_bc_type( BoundaryID bc_id, int bc_type )
  {
    _neumann_bc_map[bc_id] = bc_type;
    return;
  }

  void BCHandlingBase::set_dirichlet_bc_value( BoundaryID bc_id, Real value,
					       int component )
  {
    _dirichlet_values[bc_id](component) = value;
    return;
  }

  Real BCHandlingBase::get_dirichlet_bc_value( BoundaryID bc_id, int component ) const
  {
    return (_dirichlet_values.find(bc_id)->second)(component);
  }

  void BCHandlingBase::set_neumann_bc_value( BoundaryID bc_id, const libMesh::Point& q_in )
  {
    _q_values[bc_id] = q_in;
  }

  int BCHandlingBase::string_to_int( const std::string& bc_type_in ) const
  {
    int bc_type_out;
    if( bc_type_in == "periodic" )
      bc_type_out = PERIODIC;

    else
      {
	std::cerr << "=========================================================="  << std::endl
		  << "Error: Invalid bc_type " << bc_type_in                       << std::endl
		  << "       Physics class is " << _physics_name                   << std::endl
		  << "=========================================================="  << std::endl;
	libmesh_error();
      }
    return bc_type_out;
  }

  void BCHandlingBase::init_bc_types( const BoundaryID bc_id, 
				      const std::string& bc_id_string, 
				      const int bc_type, 
				      const GetPot& input )
  {
    switch(bc_type)
      {
      case(PERIODIC):
	{
	  /* We assume the periodic boundary pair will be listed by the master bc id.
	     Thus, if we have a periodic id and we don't see a list of bc ids with the
	     under the current id, we assume it's the slave id. We'll do consistency 
	     checks later. */
	  _num_periodic_bcs += 1;
	  int pbc_size = input.vector_variable_size("Physics/"+_physics_name+"/periodic_wall_"+bc_id_string );
	  if( pbc_size == 0 ) break;

	  // We'd better have only 2 bc ids if this is the bc list
	  if( pbc_size != 2 )
	    {
	      std::cerr << "==========================================================" 
			<< "Error: Must specify exactly two boundary condition ids " << std::endl
			<< "       for periodic boundary conditions. Found " << pbc_size << std::endl
			<< "==========================================================" << std::endl;
	      libmesh_error();
	    }
	
	  int id0 = input( "Physics/"+_physics_name+"/periodic_wall_"+bc_id_string, -1, 0 );
	  int id1 = input( "Physics/"+_physics_name+"/periodic_wall_"+bc_id_string, -1, 1 );

	  if( id0 == -1 || id1 == -1 )
	    {
	      std::cerr << "==========================================================" 
			<< "Error: Default bc id detected for periodic bc." << std::endl
			<< "       Please explicitly set periodic bc ids." << std::endl
			<< "       Detected ids " << id0 << ", " << id1 << std::endl
			<< "       for bc id = " << bc_id << std::endl
			<< "==========================================================" << std::endl;
	      libmesh_error();
	    }

	  PBCContainer pbc;
	
	  if( id0 == bc_id )
	    {
	      pbc.set_master_bcid( id0 );
	      pbc.set_slave_bcid( id1 );
	    }
	  else if( id1 == bc_id )
	    {
	      pbc.set_master_bcid( id1 );
	      pbc.set_slave_bcid( id0 );
	    }
	  else
	    {
	      // At least one of them had better be the master id
	      std::cerr << "==========================================================" 
			<< "Error: At least one of the bcs must be the master id." << std::endl
			<< "       Detected master id = " << bc_id << std::endl
			<< "       Found bc ids (" << id0 << "," << id1 << ")." << std::endl
			<< "==========================================================" << std::endl;
	      libmesh_error();
	    }

	  // Now populate offset vector
	  {
	    int offset_size = input.vector_variable_size("Physics/"+_physics_name+"/periodic_offset_"+bc_id_string );
	    if( offset_size == 0 )
	      {
		// User needs to set the offset vector for the periodic boundary
		std::cerr << "==========================================================" 
			  << "Error: Offset vector not found for bc id " << bc_id << std::endl
			  << "==========================================================" << std::endl;
		libmesh_error();
	      }
	    libMesh::RealVectorValue offset_vector;
	    for( int i = 0; i < offset_size; i++ )
	      {
		offset_vector(i) = input("Physics/"+_physics_name+"/periodic_offset_"+bc_id_string, 0.0, i );
	      }

	    pbc.set_offset_vector( offset_vector );
	  }

	  // Now stash the container object for use later in initialization
	  _periodic_bcs.push_back( pbc );

	}
	break;

      default:
	{
	  std::cerr << "==========================================================" 
		    << "Error: Invalid BC type for " << _physics_name << std::endl
		    << "       Detected BC type was " << bc_type << std::endl
		    << "==========================================================" << std::endl;
	  libmesh_error();
	}
      }
    return;
  }

  void BCHandlingBase::user_init_dirichlet_bcs( libMesh::FEMSystem* /*system*/,
						libMesh::DofMap& /*dof_map*/,
						BoundaryID /*bc_id*/,
						BCType /*bc_type*/ ) const
  {
    // Not all Physics need this so we have a do nothing default.
    return;
  }

  void BCHandlingBase::user_apply_neumann_bcs( libMesh::FEMContext& /*context*/,
					       const GRINS::CachedValues& /*cache*/,
					       const bool /*request_jacobian*/,
					       const GRINS::BoundaryID /*bc_id*/,
					       const GRINS::BCType /*bc_type*/ ) const
  {
    // Not all Physics need this so we have a do nothing default.
    return;
  }

} // namespace GRINS
