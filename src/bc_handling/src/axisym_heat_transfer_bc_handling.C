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
#include "grins/axisym_heat_transfer_bc_handling.h"

// libMesh
#include "libmesh/const_function.h"
#include "libmesh/fem_context.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  AxisymmetricHeatTransferBCHandling::AxisymmetricHeatTransferBCHandling( const std::string& physics_name,
									  const GetPot& input)
    : BCHandlingBase(physics_name)
  {
    _T_var_name = input("Physics/VariableNames/Temperature", T_var_name_default );

    std::string id_str = "Physics/"+_physics_name+"/bc_ids";
    std::string bc_str = "Physics/"+_physics_name+"/bc_types";

    this->read_bc_data( input, id_str, bc_str );

    return;
  }

  AxisymmetricHeatTransferBCHandling::~AxisymmetricHeatTransferBCHandling()
  {
    return;
  }

  int AxisymmetricHeatTransferBCHandling::string_to_int( const std::string& bc_type ) const
  {
    AHT_BC_TYPES bc_type_out;

    if( bc_type == "isothermal_wall" )
      bc_type_out = ISOTHERMAL_WALL;
  
    else if( bc_type == "adiabatic_wall" )
      bc_type_out = ADIABATIC_WALL;
  
    else if( bc_type == "prescribed_heat_flux" )
      bc_type_out = PRESCRIBED_HEAT_FLUX;
  
    else if( bc_type == "general_heat_flux" )
      bc_type_out = GENERAL_HEAT_FLUX;

    else if( bc_type == "axisymmetric" )
      {
	bc_type_out = AXISYMMETRIC;
      }
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

  void AxisymmetricHeatTransferBCHandling::init_bc_data( const libMesh::FEMSystem& system )
  {
    _T_var = system.variable_number( _T_var_name );
  }

  void AxisymmetricHeatTransferBCHandling::init_bc_types( const BoundaryID bc_id, 
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
      case(ISOTHERMAL_WALL):
	{
	  this->set_dirichlet_bc_type( bc_id, bc_type );

	  this->set_dirichlet_bc_value( bc_id, input("Physics/"+_physics_name+"/T_wall_"+bc_id_string, 0.0 ) );
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
	  std::cerr << "==========================================================" 
		    << "Error: Invalid BC type for " << _physics_name << std::endl
		    << "       Detected BC type was " << bc_type << std::endl
		    << "==========================================================" << std::endl;
	  libmesh_error();
	}
      }// End switch(bc_type)

    return;
  }

  void AxisymmetricHeatTransferBCHandling::user_apply_neumann_bcs( libMesh::FEMContext& context,
								   const GRINS::CachedValues& cache,
								   const bool request_jacobian,
								   const BoundaryID bc_id,
								   const BCType bc_type ) const
  {
    switch( bc_type )
      {
      case(AXISYMMETRIC):
	// Don't need to do anything for dT/dr = 0
	break;
	// Zero heat flux
      case(ADIABATIC_WALL):
	// Don't need to do anything: q = 0 in this case
	break;
      
	// Prescribed constant heat flux
      case(PRESCRIBED_HEAT_FLUX):
	{
	  _bound_conds.apply_neumann_axisymmetric( context, _T_var, -1.0,
						   this->get_neumann_bc_value(bc_id) );
	}
	break;
	// General heat flux from user specified function
      case(GENERAL_HEAT_FLUX):
	{
	  _bound_conds.apply_neumann_axisymmetric( context, cache, request_jacobian, _T_var, -1.0, 
						   this->get_neumann_bound_func( bc_id, _T_var ) );
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

  void AxisymmetricHeatTransferBCHandling::user_init_dirichlet_bcs( libMesh::FEMSystem* /*system*/,
								    libMesh::DofMap& dof_map,
								    BoundaryID bc_id,
								    BCType bc_type ) const
  {
    switch( bc_type )
      {
      case(ISOTHERMAL_WALL):
	{
	  std::set<BoundaryID> dbc_ids;
	  dbc_ids.insert(bc_id);
	
	  std::vector<VariableIndex> dbc_vars;
	  dbc_vars.push_back(_T_var);
	
	  ConstFunction<Number> t_func(this->get_dirichlet_bc_value(bc_id));
	
	  libMesh::DirichletBoundary t_dbc( dbc_ids, dbc_vars, &t_func );
	
	  dof_map.add_dirichlet_boundary( t_dbc );
	}
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

} // namespace GRINS
