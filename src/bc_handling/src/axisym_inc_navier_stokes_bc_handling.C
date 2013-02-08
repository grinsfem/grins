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
#include "grins/axisym_inc_navier_stokes_bc_handling.h"

//libMesh
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  AxisymmetricIncompressibleNavierStokesBCHandling::AxisymmetricIncompressibleNavierStokesBCHandling(const std::string& physics_name,
												     const GetPot& input)
    : BCHandlingBase(physics_name)
  {
    _u_r_var_name = input("Physics/VariableNames/r_velocity", u_r_var_name_default );
    _u_z_var_name = input("Physics/VariableNames/z_velocity", v_var_name_default );

    std::string id_str = "Physics/"+_physics_name+"/bc_ids";
    std::string bc_str = "Physics/"+_physics_name+"/bc_types";

    this->read_bc_data( input, id_str, bc_str );

    return;
  }

  AxisymmetricIncompressibleNavierStokesBCHandling::~AxisymmetricIncompressibleNavierStokesBCHandling()
  {
    return;
  }

  int AxisymmetricIncompressibleNavierStokesBCHandling::string_to_int( const std::string& bc_type ) const
  {
    INS_BC_TYPES bc_type_out;

    if( bc_type == "no_slip" )
      bc_type_out = NO_SLIP;

    else if( bc_type == "prescribed_vel" )
      bc_type_out = PRESCRIBED_VELOCITY;

    else if( bc_type == "general_velocity" )
      bc_type_out = GENERAL_VELOCITY;

    else if( bc_type == "axisymmetric" )
      bc_type_out = AXISYMMETRIC;

    else
      {
	std::cerr << "=========================================================="  << std::endl
		  << "Error: Invalid bc_type " << bc_type                          << std::endl
		  << "       Physics class is " << _physics_name                   << std::endl
		  << "=========================================================="  << std::endl;
	libmesh_error();
      }

    return bc_type_out;
  }
  
  void AxisymmetricIncompressibleNavierStokesBCHandling::init_bc_data( const libMesh::FEMSystem& system )
  {
    _u_r_var = system.variable_number( _u_r_var_name );
    _u_z_var = system.variable_number( _u_z_var_name );

    return;
  }

  void AxisymmetricIncompressibleNavierStokesBCHandling::init_bc_types( const BoundaryID bc_id, 
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
	
	  /* Force the user to specify 2 velocity components regardless of dimension. */
	  int n_vel_comps = input.vector_variable_size("Physics/"+_physics_name+"/bound_vel_"+bc_id_string);
	  if( n_vel_comps != 2 )
	    {
	      std::cerr << "Error: Must specify 2 velocity components when inputting"
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
	}
	break;
      case(GENERAL_VELOCITY):
	{
	  this->set_dirichlet_bc_type( bc_id, bc_type );
	}
	break;
      case(AXISYMMETRIC):
	{
	  this->set_dirichlet_bc_type( bc_id, bc_type );
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
      } // End switch(bc_type)
  
    return;
  }

  void AxisymmetricIncompressibleNavierStokesBCHandling::user_init_dirichlet_bcs( libMesh::FEMSystem* /*system*/,
										  libMesh::DofMap& dof_map,
										  BoundaryID bc_id,
										  BCType bc_type ) const
  {
    switch( bc_type )
      {
      case(NO_SLIP):
	{
	  std::set<BoundaryID> dbc_ids;
	  dbc_ids.insert(bc_id);
	
	  std::vector<VariableIndex> dbc_vars;
	  dbc_vars.push_back(_u_r_var);
	  dbc_vars.push_back(_u_z_var);
	
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
	    dbc_vars.push_back(_u_r_var);
	    ConstFunction<Number> vel_func( this->get_dirichlet_bc_value(bc_id,0) );
	  
	    libMesh::DirichletBoundary vel_dbc(dbc_ids, 
					       dbc_vars, 
					       &vel_func );
	  
	    dof_map.add_dirichlet_boundary( vel_dbc );
	    dbc_vars.clear();
	  }
	
	  {
	    dbc_vars.push_back(_u_z_var);
	    ConstFunction<Number> vel_func( this->get_dirichlet_bc_value(bc_id,1) );
	  
	    libMesh::DirichletBoundary vel_dbc(dbc_ids, 
					       dbc_vars, 
					       &vel_func );
	  
	    dof_map.add_dirichlet_boundary( vel_dbc );
	    dbc_vars.clear();
	  }
	}
	break;
      case(AXISYMMETRIC):
	{
	  std::set<BoundaryID> dbc_ids;
	  dbc_ids.insert(bc_id);
	
	  std::vector<VariableIndex> dbc_vars;
	  dbc_vars.push_back(_u_r_var);
	
	  ZeroFunction<Number> zero;
	
	  libMesh::DirichletBoundary no_slip_dbc(dbc_ids, 
						 dbc_vars, 
						 &zero );
	
	  dof_map.add_dirichlet_boundary( no_slip_dbc );
	}
	break;
      case(GENERAL_VELOCITY):
	// This case is handled in the BoundaryConditionFactory classes.
	break;
      default:
	{
	  std::cerr << "Invalid BCType " << bc_type << std::endl;
	  libmesh_error();
	}
      
      }// end switch

    return;
  }

} // namespace GRINS
