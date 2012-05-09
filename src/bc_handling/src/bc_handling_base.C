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

#include "bc_handling_base.h"

GRINS::BCHandlingBase::BCHandlingBase()
{
  return;
}

GRINS::BCHandlingBase::~BCHandlingBase()
{
  return;
}

void GRINS::BCHandlingBase::attach_neumann_bound_func( GRINS::NBCContainer& neumann_bcs )
{
  _neumann_bound_funcs = neumann_bcs;
  return;
}

void GRINS::BCHandlingBase::attach_dirichlet_bound_func( const GRINS::DBCContainer& dirichlet_bc )
{
  _dirichlet_bound_funcs.push_back( dirichlet_bc );
  return;
}

void GRINS::BCHandlingBase::read_bc_data( const GetPot& input, const std::string& id_str,
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

      this->init_bc_data( bc_id, bc_id_string, bc_type, input );
    }

  return;
}

void GRINS::BCHandlingBase::init_dirichlet_bc_func_objs( libMesh::FEMSystem* system ) const
{
  libMesh::DofMap& dof_map = system->get_dof_map();

  for( std::vector< GRINS::DBCContainer >::const_iterator 
	 it = _dirichlet_bound_funcs.begin();
       it != _dirichlet_bound_funcs.end();
       it++ )
    {
      // First, get variable names and convert them to variable id's
      std::vector<GRINS::VariableName> var_names = (*it).get_var_names();
      
      std::vector<GRINS::VariableIndex> dbc_vars;

      for( std::vector<GRINS::VariableName>::const_iterator name = var_names.begin();
	   name != var_names.end();
	   name++ )
	{
	  dbc_vars.push_back( system->variable_number( *name ) );
	}
      
      // Get bc_id's
      std::set<GRINS::BoundaryID> bc_ids = (*it).get_bc_ids();
      
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

void GRINS::BCHandlingBase::init_dirichlet_bcs( libMesh::FEMSystem* system ) const
{
  libMesh::DofMap& dof_map = system->get_dof_map();

  for( std::map< GRINS::BoundaryID,GRINS::BCType >::const_iterator it = _dirichlet_bc_map.begin();
       it != _dirichlet_bc_map.end();
       it++ )
    {
      this->user_init_dirichlet_bcs( system, dof_map, it->first, it->second );
    }

  return;
}

void GRINS::BCHandlingBase::set_dirichlet_bc_type( GRINS::BoundaryID bc_id, int bc_type )
{
  _dirichlet_bc_map[bc_id] = bc_type;
  return;
}

void GRINS::BCHandlingBase::set_neumann_bc_type( GRINS::BoundaryID bc_id, int bc_type )
{
  _neumann_bc_map[bc_id] = bc_type;
  return;
}

void GRINS::BCHandlingBase::set_dirichlet_bc_value( GRINS::BoundaryID bc_id, Real value,
						    int component )
{
  _dirichlet_values[bc_id](component) = value;
  return;
}

Real GRINS::BCHandlingBase::get_dirichlet_bc_value( GRINS::BoundaryID bc_id, int component ) const
{
  return (_dirichlet_values.find(bc_id)->second)(component);
}

void GRINS::BCHandlingBase::set_neumann_bc_value( GRINS::BoundaryID bc_id, const libMesh::Point& q_in )
{
  _q_values[bc_id] = q_in;
}

int GRINS::BCHandlingBase::string_to_int( const std::string& bc_type_in ) const
{
  // Default to negative value to help catch forgetting to overload this when
  // necessary.
  return -1;
}

void GRINS::BCHandlingBase::init_bc_data( const GRINS::BoundaryID bc_id, 
					  const std::string& bc_id_string, 
					  const int bc_type, 
					  const GetPot& input )
{
  // Not all Physics need this so we have a do nothing default.
  return;
}

void GRINS::BCHandlingBase::user_init_dirichlet_bcs( libMesh::FEMSystem* system, libMesh::DofMap& dof_map,
						     GRINS::BoundaryID bc_id, GRINS::BCType bc_type ) const
{
  // Not all Physics need this so we have a do nothing default.
  return;
}
