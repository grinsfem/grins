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

#include "low_mach_navier_stokes_bc_handling.h"

namespace GRINS
{

  LowMachNavierStokesBCHandling::LowMachNavierStokesBCHandling(const std::string& physics_name,
							       const GetPot& input)
    : BCHandlingBase(physics_name)
  {
    _u_var_name = input("Physics/VariableNames/u_velocity", u_var_name_default );
    _v_var_name = input("Physics/VariableNames/v_velocity", v_var_name_default );
    _w_var_name = input("Physics/VariableNames/w_velocity", w_var_name_default );
    _T_var_name = input("Physics/VariableNames/Temperature", T_var_name_default );

    std::string id_str = "Physics/"+_physics_name+"/vel_bc_ids";
    std::string bc_str = "Physics/"+_physics_name+"/vel_bc_types";

    this->read_bc_data( input, id_str, bc_str );

    id_str = "Physics/"+_physics_name+"/temp_bc_ids";
    bc_str = "Physics/"+_physics_name+"/temp_bc_types";

    this->read_bc_data( input, id_str, bc_str );

    return;
  }

  LowMachNavierStokesBCHandling::~LowMachNavierStokesBCHandling()
  {
    return;
  }

  int LowMachNavierStokesBCHandling::string_to_int( const std::string& bc_type ) const
  {
    int bc_type_out;

    if( bc_type == "no_slip" )
      bc_type_out = NO_SLIP;

    else if( bc_type == "prescribed_vel" )
      bc_type_out = PRESCRIBED_VELOCITY;

    else if( bc_type == "inflow" )
      bc_type_out = INFLOW;

    else if( bc_type == "isothermal" )
      bc_type_out = ISOTHERMAL_WALL;
  
    else if( bc_type == "adiabatic" )
      bc_type_out = ADIABATIC_WALL;
  
    else if( bc_type == "prescribed_heat_flux" )
      bc_type_out = PRESCRIBED_HEAT_FLUX;
  
    else if( bc_type == "general_heat_flux" )
      bc_type_out = GENERAL_HEAT_FLUX;

    else
      {
	// Call base class to detect any physics-common boundary conditions
	bc_type_out = BCHandlingBase::string_to_int( bc_type );
      }

    return bc_type_out;
  }

  void LowMachNavierStokesBCHandling::init_bc_data( const BoundaryID bc_id, 
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
      case(ISOTHERMAL_WALL):
	{
	  this->set_temp_bc_type( bc_id, bc_type );

	  this->set_temp_bc_value( bc_id, input("Physics/"+_physics_name+"/T_wall_"+bc_id_string, 0.0 ) );
	}
	break;
      
      case(ADIABATIC_WALL):
	{
	  this->set_neumann_bc_type( bc_id, bc_type );
	}
	break;
      case(PRESCRIBED_HEAT_FLUX):
	{
	  this->set_neumann_bc_type( bc_id, bc_type );
	
	  libMesh::Point q_in;
	
	  int num_q_components = input.vector_variable_size("Physics/"+_physics_name+"/q_wall_"+bc_id_string);
	
	  for( int i = 0; i < num_q_components; i++ )
	    {
	      q_in(i) = input("Physics/"+_physics_name+"/q_wall_"+bc_id_string, 0.0, i );
	    }

	  this->set_neumann_bc_value( bc_id, q_in );
	}
	break;
      case(GENERAL_HEAT_FLUX):
	{
	  this->set_neumann_bc_type( bc_id, bc_type );
	}
	break;
      default:
	{
	  // Call base class to detect any physics-common boundary conditions
	  BCHandlingBase::init_bc_data( bc_id, bc_id_string, bc_type, input );
	}
      } // End switch(bc_type)
  
    return;
  }

  void LowMachNavierStokesBCHandling::user_init_dirichlet_bcs( libMesh::FEMSystem* system,
							       libMesh::DofMap& dof_map,
							       BoundaryID bc_id,
							       BCType bc_type ) const
  {
    int dim = system->get_mesh().mesh_dimension();

    VariableIndex T_var = system->variable_number( _T_var_name );
    VariableIndex u_var = system->variable_number( _u_var_name );
    VariableIndex v_var = system->variable_number( _v_var_name );
    VariableIndex w_var;
    if( dim == 3 )
      w_var = system->variable_number( _w_var_name );

    switch( bc_type )
      {
      case(NO_SLIP):
	{
	  std::set<BoundaryID> dbc_ids;
	  dbc_ids.insert(bc_id);
	
	  std::vector<VariableIndex> dbc_vars;
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
	  std::set<BoundaryID> dbc_ids;
	  dbc_ids.insert(bc_id);
	
	  std::vector<VariableIndex> dbc_vars;
	
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
      case(ISOTHERMAL_WALL):
	{
	  std::set<BoundaryID> dbc_ids;
	  dbc_ids.insert(bc_id);
	
	  std::vector<VariableIndex> dbc_vars;
	  dbc_vars.push_back(T_var);
	
	  ConstFunction<Number> t_func(this->get_temp_bc_value(bc_id));
	
	  libMesh::DirichletBoundary t_dbc( dbc_ids, dbc_vars, &t_func );
	
	  dof_map.add_dirichlet_boundary( t_dbc );
	}
	break;
      default:
	{
	  std::cerr << "Invalid BCType " << bc_type << std::endl;
	  libmesh_error();
	}
      
      }// end switch

    return;
  }

  void LowMachNavierStokesBCHandling::set_temp_bc_type( BoundaryID bc_id, int bc_type )
  {
    _temp_bc_map[bc_id] = bc_type;
    return;
  }

  void LowMachNavierStokesBCHandling::set_temp_bc_value( BoundaryID bc_id, Real value )
  {
    _T_values[bc_id] = value;
    return;
  }
  Real LowMachNavierStokesBCHandling::get_temp_bc_value( BoundaryID bc_id ) const
  {
    return _T_values.find(bc_id)->second;
  }

  void LowMachNavierStokesBCHandling::init_dirichlet_bcs( libMesh::FEMSystem* system ) const
  {
    libMesh::DofMap& dof_map = system->get_dof_map();

    for( std::map< BoundaryID,BCType >::const_iterator it = _dirichlet_bc_map.begin();
	 it != _dirichlet_bc_map.end();
	 it++ )
      {
	this->user_init_dirichlet_bcs( system, dof_map, it->first, it->second );
      }

    for( std::map< BoundaryID,BCType >::const_iterator it = _temp_bc_map.begin();
	 it != _temp_bc_map.end();
	 it++ )
      {
	this->user_init_dirichlet_bcs( system, dof_map, it->first, it->second );
      }

    return;
  }

} // namespace GRINS
