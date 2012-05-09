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

#include "inc_navier_stokes_bc_handling.h"

GRINS::IncompressibleNavierStokesBCHandling::IncompressibleNavierStokesBCHandling(std::string& physics_name,
										  const GetPot& input)
  : BCHandlingBase(),
    _physics_name(physics_name)
{
  _u_var_name = input("Physics/VariableNames/Temperature", GRINS::u_var_name_default );
  _v_var_name = input("Physics/VariableNames/Temperature", GRINS::v_var_name_default );
  _w_var_name = input("Physics/VariableNames/Temperature", GRINS::w_var_name_default );

  std::string id_str = "Physics/"+_physics_name+"/bc_ids";
  std::string bc_str = "Physics/"+_physics_name+"/bc_types";

  this->read_bc_data( input, id_str, bc_str );

  return;
}

GRINS::IncompressibleNavierStokesBCHandling::~IncompressibleNavierStokesBCHandling()
{
  return;
}

int GRINS::IncompressibleNavierStokesBCHandling::string_to_int( const std::string& bc_type ) const
{
  INS_BC_TYPES bc_type_out;

  if( bc_type == "no_slip" )
    bc_type_out = NO_SLIP;

  else if( bc_type == "prescribed_vel" )
    bc_type_out = PRESCRIBED_VELOCITY;

  else if( bc_type == "inflow" )
    bc_type_out = INFLOW;

  else
    {
      std::cerr << "Error: Invalid bc_type " << bc_type << std::endl;
      libmesh_error();
    }

  return bc_type_out;
}

void GRINS::IncompressibleNavierStokesBCHandling::init_bc_data( const GRINS::BoundaryID bc_id, 
								const std::string& bc_id_string, 
								const int bc_type, 
								const GetPot& input )
{
  switch(bc_type)
    {
    case(NO_SLIP):
      {
	this->set_dirichlet_bc_type( bc_id, bc_type );
      }
      break;
    case(PRESCRIBED_VELOCITY):
      {
	this->set_dirichlet_bc_type( bc_id, bc_type );
	
	/* Force the user to specify 3 velocity components regardless of dimension.
	   This should make it easier to keep things correct if we want to have 
	   2D flow not be in the x-y plane. */
	int n_vel_comps = input.vector_variable_size("Physics/"+_physics_name+"/bound_vel_"+bc_id_string);
	if( n_vel_comps != 3 )
	  {
	    std::cerr << "Error: Must specify 3 velocity components when inputting"
		      << std::endl
		      << "       prescribed velocities. Found " << n_vel_comps
		      << " velocity components."
		      << std::endl;
	    libmesh_error();
	  }
	
	this->set_dirichlet_bc_value( bc_id, 
				      input("Physics/"+_physics_name+"/bound_vel_"+bc_id_string, 0.0, 0 ),
				      0 );

	this->set_dirichlet_bc_value( bc_id, 
				      input("Physics/"+_physics_name+"/bound_vel_"+bc_id_string, 0.0, 1 ),
				      1 );

	this->set_dirichlet_bc_value( bc_id, 
				      input("Physics/"+_physics_name+"/bound_vel_"+bc_id_string, 0.0, 2 ),
				      2 );
      }
      break;
    case(INFLOW):
      {
	this->set_dirichlet_bc_type( bc_id, bc_type );
      }
      break;
    default:
      {
	std::cerr << "Error: Invalid Dirichlet BC type for " << _physics_name
		  << std::endl;
	libmesh_error();
      }
    } // End switch(bc_type)
  
  return;
}

void GRINS::IncompressibleNavierStokesBCHandling::user_init_dirichlet_bcs( libMesh::FEMSystem* system,
									   libMesh::DofMap& dof_map,
									   GRINS::BoundaryID bc_id,
									   GRINS::BCType bc_type ) const
{
  int dim = system->get_mesh().mesh_dimension();

  GRINS::VariableIndex u_var = system->variable_number( _u_var_name );
  GRINS::VariableIndex v_var = system->variable_number( _v_var_name );
  GRINS::VariableIndex w_var;
  if( dim == 3 )
    w_var = system->variable_number( _w_var_name );

  switch( bc_type )
    {
    case(NO_SLIP):
      {
	std::set<GRINS::BoundaryID> dbc_ids;
	dbc_ids.insert(bc_id);
	
	std::vector<GRINS::VariableIndex> dbc_vars;
	dbc_vars.push_back(u_var);
	dbc_vars.push_back(v_var);
	if(dim == 3)
	  dbc_vars.push_back(w_var);
	
	ZeroFunction<Number> zero;
	
	libMesh::DirichletBoundary no_slip_dbc(dbc_ids, 
					       dbc_vars, 
					       &zero );
	
	dof_map.add_dirichlet_boundary( no_slip_dbc );
      }
      break;
    case(PRESCRIBED_VELOCITY):
      {
	std::set<GRINS::BoundaryID> dbc_ids;
	dbc_ids.insert(bc_id);
	
	std::vector<GRINS::VariableIndex> dbc_vars;
	
	// This is inefficient, but it shouldn't matter because
	// everything gets cached on the libMesh side so it should
	// only affect performance at startup.
	{
	  dbc_vars.push_back(u_var);
	  ConstFunction<Number> vel_func( this->get_dirichlet_bc_value(bc_id,0) );
	  
	  libMesh::DirichletBoundary vel_dbc(dbc_ids, 
					     dbc_vars, 
					     &vel_func );
	  
	  dof_map.add_dirichlet_boundary( vel_dbc );
	  dbc_vars.clear();
	}
	
	{
	  dbc_vars.push_back(v_var);
	  ConstFunction<Number> vel_func( this->get_dirichlet_bc_value(bc_id,1) );
	  
	  libMesh::DirichletBoundary vel_dbc(dbc_ids, 
					     dbc_vars, 
					     &vel_func );
	  
	  dof_map.add_dirichlet_boundary( vel_dbc );
	  dbc_vars.clear();
	}
	if( dim == 3 )
	  {
	    dbc_vars.push_back(w_var);
	    ConstFunction<Number> vel_func( this->get_dirichlet_bc_value(bc_id,2) );
	    
	    libMesh::DirichletBoundary vel_dbc(dbc_ids, 
					       dbc_vars, 
					       &vel_func );
	    
	    dof_map.add_dirichlet_boundary( vel_dbc );
	  }  
      }
      break;
    case(INFLOW):
      // This case is handled in the BoundaryConditionFactory classes.
      break;
    default:
      {
	std::cerr << "Invalid GRINS::BCType " << bc_type << std::endl;
	libmesh_error();
      }
      
    }// end switch

  return;
}
