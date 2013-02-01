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

#include "grins/axisym_electrostatics_bc_handling.h"

namespace GRINS
{

  AxisymmetricElectrostaticsBCHandling::AxisymmetricElectrostaticsBCHandling( const std::string& physics_name,
									      const GetPot& input)
    : BCHandlingBase(physics_name)
  {
    _V_var_name = input("Physics/VariableNames/ElectricPotential", GRINS::V_var_name_default );
    
    std::string id_str = "Physics/"+_physics_name+"/bc_ids";
    std::string bc_str = "Physics/"+_physics_name+"/bc_types";
    
    this->read_bc_data( input, id_str, bc_str );
    
    return;
  }

  GRINS::AxisymmetricElectrostaticsBCHandling::~AxisymmetricElectrostaticsBCHandling()
  {
    return;
  }

  int GRINS::AxisymmetricElectrostaticsBCHandling::string_to_int( const std::string& bc_type ) const
  {
    AE_BC_TYPES bc_type_out;
    
    if( bc_type == "insulating_wall" )
      bc_type_out = INSULATING_WALL;
    
    else if( bc_type == "prescribed_current" )
      bc_type_out = PRESCRIBED_CURRENT;
    
    else if( bc_type == "general_current" )
      bc_type_out = GENERAL_CURRENT;
    
    else if( bc_type == "prescribed_voltage" )
      bc_type_out = PRESCRIBED_VOLTAGE;

    else if( bc_type == "general_voltage" )
      bc_type_out = GENERAL_VOLTAGE;
    
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
  
  void GRINS::AxisymmetricElectrostaticsBCHandling::init_bc_data( const GRINS::BoundaryID bc_id, 
								  const std::string& bc_id_string, 
								  const int bc_type, 
								  const GetPot& input )
  {
  switch(bc_type)
    {
    case(AXISYMMETRIC):
      {
	this->set_neumann_bc_type( bc_id, bc_type );
      }
      break;
    case(INSULATING_WALL):
      // Do nothing BC
      break;
      
    case(PRESCRIBED_CURRENT):
      {
	this->set_neumann_bc_type( bc_id, bc_type );

	int num_E_components = input.vector_variable_size("Physics/"+_physics_name+"/E_wall_"+bc_id_string);
	
	libMesh::Point E_in;

	for( int i = 0; i < num_E_components; i++ )
	  {
	    E_in(i) = input("Physics/"+_physics_name+"/E_wall_"+bc_id_string, 0.0, i );
	  }

	this->set_neumann_bc_value( bc_id, E_in );
      }
      break;
      
    case(GENERAL_CURRENT):
      {
	this->set_neumann_bc_type( bc_id, bc_type );
      }
      break;

    case(PRESCRIBED_VOLTAGE):
      {
	this->set_dirichlet_bc_type( bc_id, bc_type );

	this->set_dirichlet_bc_value( bc_id, input("Physics/"+_physics_name+"/V_wall_"+bc_id_string, 0.0 ) );
      }
      break;

    case(GENERAL_VOLTAGE):
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
    }// End switch(bc_type)

  return;
}

  void GRINS::AxisymmetricElectrostaticsBCHandling::user_apply_neumann_bcs( libMesh::FEMContext& context,
									    GRINS::VariableIndex var,
									    bool request_jacobian,
									    GRINS::BoundaryID bc_id,
									    GRINS::BCType bc_type ) const
  {
    switch( bc_type )
      {
      case(AXISYMMETRIC):
	// Don't need to do anything for dV/dr = 0
	break;

      case(INSULATING_WALL):
	// Do nothing BC
	break;
	
      case(PRESCRIBED_CURRENT):
	{
	  _bound_conds.apply_neumann_axisymmetric( context, var, -1.0,
						   this->get_neumann_bc_value(bc_id) );
	}
	break;
	
      case(GENERAL_CURRENT):
	{
	  _bound_conds.apply_neumann_axisymmetric( context, request_jacobian, var, -1.0, 
						   this->get_neumann_bound_func( bc_id, var ) );
	}
	break;
	
      default:
	{
	  std::cerr << "Error: Invalid Neumann BC type for " << _physics_name
		    << std::endl;
	  libmesh_error();
	}
      }
    return;
  }
  
  void GRINS::AxisymmetricElectrostaticsBCHandling::user_init_dirichlet_bcs( libMesh::FEMSystem* system,
									     libMesh::DofMap& dof_map,
									     GRINS::BoundaryID bc_id,
									     GRINS::BCType bc_type ) const
  {
    GRINS::VariableIndex V_var = system->variable_number( _V_var_name );
    
    switch( bc_type )
      {

      case(PRESCRIBED_VOLTAGE):
	{
	  std::set<GRINS::BoundaryID> dbc_ids;
	  dbc_ids.insert(bc_id);
	
	  std::vector<GRINS::VariableIndex> dbc_vars;
	  dbc_vars.push_back(V_var);

	  ConstFunction<Number> voltage_func( this->get_dirichlet_bc_value(bc_id) );
	  
	  libMesh::DirichletBoundary voltage_dbc(dbc_ids, 
						 dbc_vars, 
						 &voltage_func );
	  
	  dof_map.add_dirichlet_boundary( voltage_dbc );
	}

      case(GENERAL_CURRENT):
	// This case is handled in the BoundaryConditionFactory classes.
	break;

      default:
	{
	  std::cerr << "Error: Invalid Dirichlet BC type for " << _physics_name
		    << std::endl;
	  libmesh_error();
	}
      }// end switch
    
    return;
  }


} //namespace GRINS
