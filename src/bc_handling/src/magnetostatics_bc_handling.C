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
#include "grins/magnetostatics_bc_handling.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/point.h"
#include "libmesh/const_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/fem_system.h"
#include "libmesh/fem_context.h"

namespace GRINS
{

  MagnetostaticsBCHandling::MagnetostaticsBCHandling( const std::string& physics_name,
						      const GetPot& input)
    : BCHandlingBase(physics_name)
  {
    _A_var_name = input("Physics/VariableNames/MagneticPotential", GRINS::A_var_name_default );
    
    std::string id_str = "Physics/"+_physics_name+"/bc_ids";
    std::string bc_str = "Physics/"+_physics_name+"/bc_types";
    std::string var_str = "Physics/"+_physics_name+"/bc_variables";
    std::string val_str = "Physics/"+_physics_name+"/bc_values";

    this->read_bc_data( input, id_str, bc_str, var_str, val_str );
    
    return;
  }

  MagnetostaticsBCHandling::~MagnetostaticsBCHandling()
  {
    return;
  }

  void MagnetostaticsBCHandling::init_bc_data( const libMesh::FEMSystem& system )
  {
    _A_var = system.variable_number(_A_var_name);

    return;
  }

  int MagnetostaticsBCHandling::string_to_int( const std::string& bc_type ) const
  {
    M_BC_TYPES bc_type_out;
    
    if( bc_type == "unbounded" )
      bc_type_out = UNBOUNDED;

    else if( bc_type == "prescribed_magnetic_flux" )
      bc_type_out = PRESCRIBED_MAGNETIC_FLUX;

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
  
  void MagnetostaticsBCHandling::init_bc_types( const GRINS::BoundaryID bc_id, 
                                                const std::string& bc_id_string, 
                                                const int bc_type, 
                                                const std::string& bc_vars, 
                                                const std::string& bc_value, 
                                                const GetPot& input )
  {
    switch(bc_type)
      {
      case(UNBOUNDED):
	// Do nothing BC
	break;
  
      case(PRESCRIBED_MAGNETIC_FLUX):
	{
	  this->set_neumann_bc_type( bc_id, bc_type );

	  libMesh::Point B_in;
	
	  int num_B_components = input.vector_variable_size("Physics/"+_physics_name+"/B_wall_"+bc_id_string);
	
	  for( int i = 0; i < num_B_components; i++ )
	    {
	      B_in(i) = input("Physics/"+_physics_name+"/B_wall_"+bc_id_string, 0.0, i );
	    }

	  this->set_neumann_bc_value( bc_id, B_in );
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

  void MagnetostaticsBCHandling::user_apply_neumann_bcs( AssemblyContext& context,
							 const GRINS::CachedValues& cache,
							 bool request_jacobian,
							 BoundaryID bc_id,
							 BCType bc_type ) const
  {
    switch( bc_type )
      {
      case(PRESCRIBED_MAGNETIC_FLUX):
	{
	  _bound_conds.apply_neumann_cross( context, _A_var, 1.0,
					    this->get_neumann_bc_value(bc_id) );
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
  
  void MagnetostaticsBCHandling::user_init_dirichlet_bcs( libMesh::FEMSystem* system,
							  libMesh::DofMap& dof_map,
							  BoundaryID bc_id,
							  BCType bc_type ) const
  {
    VariableIndex A_var = system->variable_number( _A_var_name );
    
    switch( bc_type )
      {
	
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
