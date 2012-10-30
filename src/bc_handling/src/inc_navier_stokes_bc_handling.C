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

GRINS::IncompressibleNavierStokesBCHandling::IncompressibleNavierStokesBCHandling(const std::string& physics_name,
										  const GetPot& input)
  : GRINS::BCHandlingBase(physics_name)
{
  _u_var_name = input("Physics/VariableNames/u_velocity", GRINS::u_var_name_default );
  _v_var_name = input("Physics/VariableNames/v_velocity", GRINS::v_var_name_default );
  _w_var_name = input("Physics/VariableNames/w_velocity", GRINS::w_var_name_default );

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
  int bc_type_out;

  if( bc_type == "no_slip" )
    bc_type_out = NO_SLIP;

  else if( bc_type == "prescribed_vel" )
    bc_type_out = PRESCRIBED_VELOCITY;

  else if( bc_type == "parabolic_profile" )
    bc_type_out = PARABOLIC_PROFILE;

  else if( bc_type == "inflow" )
    bc_type_out = INFLOW;

  else
    {
      // Call base class to detect any physics-common boundary conditions
      bc_type_out = GRINS::BCHandlingBase::string_to_int( bc_type );
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
    case(PARABOLIC_PROFILE):
      {
	this->set_dirichlet_bc_type( bc_id, bc_type );
	
	// Make sure all 6 components are there
	if( input.vector_variable_size("Physics/"+_physics_name+"/parabolic_profile_coeffs_"+bc_id_string) != 6 )
	  {
	    std::cerr << "Error: Must specify 6 components when inputting"
		      << std::endl
		      << "       coefficients for a parabolic profile. Found " 
		      << input.vector_variable_size("Physics/"+_physics_name+"/parabolic_profile_"+bc_id_string)
		      << " components."
		      << std::endl;
	    libmesh_error();
	  }

	Real a = input( "Physics/"+_physics_name+"/parabolic_profile_coeffs_"+bc_id_string, 0.0, 0 );
	Real b = input( "Physics/"+_physics_name+"/parabolic_profile_coeffs_"+bc_id_string, 0.0, 1 );
	Real c = input( "Physics/"+_physics_name+"/parabolic_profile_coeffs_"+bc_id_string, 0.0, 2 );
	Real d = input( "Physics/"+_physics_name+"/parabolic_profile_coeffs_"+bc_id_string, 0.0, 3 );
	Real e = input( "Physics/"+_physics_name+"/parabolic_profile_coeffs_"+bc_id_string, 0.0, 4 );
	Real f = input( "Physics/"+_physics_name+"/parabolic_profile_coeffs_"+bc_id_string, 0.0, 5 );

	std::string var = input( "Physics/"+_physics_name+"/parabolic_profile_var_"+bc_id_string, "DIE!" );
	
	GRINS::DBCContainer cont;
	cont.add_var_name( var );
	cont.add_bc_id( bc_id );

	std::tr1::shared_ptr<libMesh::FunctionBase<Number> > func( new GRINS::ParabolicProfile(a,b,c,d,e,f) );
	cont.set_func( func );
	this->attach_dirichlet_bound_func( cont );
	
	// Set specified components of Dirichlet data to zero
	std::string fix_var = input( "Physics/"+_physics_name+"/parabolic_profile_fix_"+bc_id_string, "DIE!" );

	GRINS::DBCContainer cont_fix;
	cont_fix.add_var_name( fix_var );
	cont_fix.add_bc_id( bc_id );

	std::tr1::shared_ptr<libMesh::FunctionBase<Number> > func_fix( new ZeroFunction<Number>() );
	cont_fix.set_func( func_fix );
	this->attach_dirichlet_bound_func( cont_fix );
      }
      break;

    case(INFLOW):
      {
	this->set_dirichlet_bc_type( bc_id, bc_type );
      }
      break;

    default:
      {
	// Call base class to detect any physics-common boundary conditions
	GRINS::BCHandlingBase::init_bc_data( bc_id, bc_id_string, bc_type, input );
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
    case(PARABOLIC_PROFILE):
      // This case is handled init_dirichlet_bc_func_objs
      break;
    case(INFLOW):
      // This case is handled in the init_dirichlet_bc_func_objs
      break;
    default:
      {
	std::cerr << "Invalid GRINS::BCType " << bc_type << std::endl;
	libmesh_error();
      }
      
    }// end switch

  return;
}
