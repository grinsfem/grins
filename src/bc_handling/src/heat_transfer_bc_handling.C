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

#include "heat_transfer_bc_handling.h"

GRINS::HeatTransferBCHandling::HeatTransferBCHandling(std::string& physics_name)
  : BCHandlingBase(),
    _physics_name(physics_name)
{
  return;
}

GRINS::HeatTransferBCHandling::~HeatTransferBCHandling()
{
  return;
}

int GRINS::HeatTransferBCHandling::string_to_int( const std::string& bc_type ) const
{
  HT_BC_TYPES bc_type_out;

  if( bc_type == "isothermal_wall" )
    bc_type_out = ISOTHERMAL_WALL;
  
  else if( bc_type == "adiabatic_wall" )
    bc_type_out = ADIABATIC_WALL;
  
  else if( bc_type == "prescribed_heat_flux" )
    bc_type_out = PRESCRIBED_HEAT_FLUX;
  
  else if( bc_type == "general_heat_flux" )
    bc_type_out = GENERAL_HEAT_FLUX;

  else
    {
      std::cerr << "Error: Invalid bc_type " << bc_type << std::endl;
      libmesh_error();
    }

  return bc_type_out;
}

void GRINS::HeatTransferBCHandling::init_bc_data( const GRINS::BoundaryID bc_id, 
						  const std::string& bc_id_string, 
						  const int bc_type, 
						  const GetPot& input )
{
  switch(bc_type)
    {
    case(ISOTHERMAL_WALL):
      {
	_dirichlet_bc_map[bc_id] = bc_type;
	
	_dirichlet_values[bc_id] = 
	  input("Physics/"+_physics_name+"/T_wall_"+bc_id_string, 0.0 );
      }
      break;
      
    case(ADIABATIC_WALL):
      {
	_neumann_bc_map[bc_id] = bc_type;
      }
      break;
      
    case(PRESCRIBED_HEAT_FLUX):
      {
	_neumann_bc_map[bc_id] = bc_type;
	
	libMesh::Point q_in;
	
	int num_q_components = input.vector_variable_size("Physics/"+_physics_name+"/q_wall_"+bc_id_string);
	
	for( int i = 0; i < num_q_components; i++ )
	  {
	    q_in(i) = input("Physics/"+_physics_name+"/q_wall_"+bc_id_string, 0.0, i );
	  }
	_q_values[bc_id] = q_in;
      }
      break;
    case(GENERAL_HEAT_FLUX):
      {
	_neumann_bc_map[bc_id] = bc_type;
      }
      break;
    default:
      {
	std::cerr << "Error: Invalid Dirichlet BC type for " << _physics_name
		  << std::endl;
	libmesh_error();
      }
      
    }// End switch(bc_type)

  return;
}

void GRINS::HeatTransferBCHandling::init_dirichlet_bcs( libMesh::DofMap& dof_map ) const
{
  for( std::map< GRINS::BoundaryID,GRINS::BCType >::const_iterator it = _dirichlet_bc_map.begin();
       it != _dirichlet_bc_map.end();
       it++ )
    {
      switch( it->second )
	{
	case(ISOTHERMAL_WALL):
	  {
	    std::set<GRINS::BoundaryID> dbc_ids;
	    dbc_ids.insert(it->first);

	    std::vector<GRINS::VariableIndex> dbc_vars;
	    //	    dbc_vars.push_back(_T_var);

	    //Real value = _dirichlet_values.find(it->first);

	    //ConstFunction<Number> t_func( _dirichlet_values[it->first] );

	    //libMesh::DirichletBoundary t_dbc( dbc_ids, dbc_vars, &t_func );

	    //dof_map.add_dirichlet_boundary( t_dbc );
	  }
	  break;
	default:
	  {
	    std::cerr << "Error: Invalid Dirichlet BC type for " << _physics_name
		      << std::endl;
	    libmesh_error();
	  }
  	}// end switch
    } //end for

  return;
}
