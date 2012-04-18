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

GRINS::Physics::Physics( const std::string& physics_name )
  : _physics_name( physics_name )
{
  return;
}

GRINS::Physics::~Physics()
{
  return;
}

void GRINS::Physics::read_input_options( const GetPot& input )
{
  return;
}

void GRINS::Physics::set_time_evolving_vars( libMesh::FEMSystem* system )
{
  return;
}

void GRINS::Physics::read_bc_data( const GetPot& input )
{
  int num_ids = input.vector_variable_size("Physics/"+_physics_name+"/bc_ids");
  int num_bcs = input.vector_variable_size("Physics/"+_physics_name+"/bc_types");

  if( num_ids != num_bcs )
    {
      std::cerr << "Error: Must specify equal number of boundary ids and boundary conditions"
		<< std::endl;
      libmesh_error();
    }

  for( int i = 0; i < num_ids; i++ )
    {
      int bc_id = input("Physics/"+_physics_name+"/bc_ids", -1, i );
      std::string bc_type_in = input("Physics/"+_physics_name+"/bc_types", "NULL", i );

      int bc_type = this->string_to_int( bc_type_in );

      std::stringstream ss;
      ss << bc_id;
      std::string bc_id_string = ss.str();

      this->init_bc_data( bc_id, bc_id_string, bc_type, input );
    }

  return;
}

int GRINS::Physics::string_to_int( const std::string& bc_type_in )
{
  // Default to negative value to help catch forgetting to overload this when
  // necessary.
  return -1;
}

void GRINS::Physics::init_bc_data( const GRINS::BoundaryID bc_id, 
				   const std::string& bc_id_string, 
				   const int bc_type, 
				   const GetPot& input )
{
  // Not all Physics need this so we have a do nothing default.
  return;
}

void GRINS::Physics::init_dirichlet_bcs( libMesh::DofMap& dof_map )
{
  // Not all Physics need this so we have a do nothing default.
  return;
}

void GRINS::Physics::init_user_dirichlet_bcs( libMesh::FEMSystem* system )
{
  libMesh::DofMap& dof_map = system->get_dof_map();

  for( std::vector< GRINS::DBCContainer >::iterator 
	 it = _dirichlet_bound_funcs.begin();
       it != _dirichlet_bound_funcs.end();
       it++ )
    {
      std::vector<GRINS::VariableName> var_names = (*it).get_var_names();
      
      std::vector<GRINS::VariableIndex> dbc_vars;

      for( std::vector<GRINS::VariableName>::const_iterator name = var_names.begin();
	   name != var_names.end();
	   name++ )
	{
	  dbc_vars.push_back( system->variable_number( *name ) );
	}
      
      std::set<GRINS::BoundaryID> bc_ids = (*it).get_bc_ids();
      
      std::tr1::shared_ptr<libMesh::FunctionBase<Number> > func = (*it).get_func();

      libMesh::DirichletBoundary dbc( bc_ids, dbc_vars, &*func );
      
      dof_map.add_dirichlet_boundary( dbc );
    }

  return;
}

void GRINS::Physics::attach_neumann_bound_func( GRINS::NBCContainer& neumann_bcs )
{
  _neumann_bound_funcs = neumann_bcs;

  return;
}

void GRINS::Physics::attach_dirichlet_bound_func( const GRINS::DBCContainer& dirichlet_bc )
{
  _dirichlet_bound_funcs.push_back( dirichlet_bc );
  return;
}

#ifdef USE_GRVY_TIMERS
void GRINS::Physics::attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer )
{
  _timer = grvy_timer;
  return;
}
#endif
