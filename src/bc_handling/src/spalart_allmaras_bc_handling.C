//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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


// This class
#include "grins/spalart_allmaras_bc_handling.h"

// libMesh
#include "libmesh/const_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  SpalartAllmarasBCHandling::SpalartAllmarasBCHandling(const std::string& physics_name,
						 const GetPot& input)
    : BCHandlingBase(physics_name),
      _turb_vars(input)
  {
    std::string id_str = "Physics/"+_physics_name+"/bc_ids";
    std::string bc_str = "Physics/"+_physics_name+"/bc_types";
    std::string var_str = "Physics/"+_physics_name+"/bc_variables";
    std::string val_str = "Physics/"+_physics_name+"/bc_values";

    this->read_bc_data( input, id_str, bc_str, var_str, val_str );
    
    return;
  }

  SpalartAllmarasBCHandling::~SpalartAllmarasBCHandling()
  {
    return;
  }

  int SpalartAllmarasBCHandling::string_to_int( const std::string& bc_type ) const
  {
    int bc_type_out;

    if( bc_type == "general_velocity" )
      {
	bc_type_out = GENERAL_VELOCITY;
      }
    else
      {
	// Call base class to detect any physics-common boundary conditions
	bc_type_out = BCHandlingBase::string_to_int( bc_type );
      }

    return bc_type_out;
  }

  void SpalartAllmarasBCHandling::init_bc_data( const libMesh::FEMSystem& system )
  {
    _turb_vars.init(const_cast<libMesh::FEMSystem*>(&system));

    return;
  }
  
  void SpalartAllmarasBCHandling::init_bc_types( const BoundaryID bc_id, 
					      const std::string& bc_id_string, 
					      const int bc_type, 
					      const std::string& bc_vars, 
					      const std::string& bc_value, 
					      const GetPot& input )
  { 
    switch(bc_type)
      {
      case(GENERAL_VELOCITY):
	{
	  this->set_dirichlet_bc_type( bc_id, bc_type);
	}
	break;
	
      default:
	{
	  // Call base class to detect any physics-common boundary conditions
	  BCHandlingBase::init_bc_types( bc_id, bc_id_string, bc_type,
                                         bc_vars, bc_value, input );	
	}
      } // End switch(bc_type)
    return;
  }
  
   void IncompressibleNavierStokesBCHandling::user_init_dirichlet_bcs( libMesh::FEMSystem* system,
								      libMesh::DofMap& dof_map,
								      BoundaryID bc_id,
								      BCType bc_type ) const
  {
    switch( bc_type )
      {
	case(GENERAL_VELOCITY):
	// This case is handled in the init_dirichlet_bc_func_objs
	break;
	
      default:
	{
	  std::cerr << "Invalid BCType " << bc_type << std::endl;
	  libmesh_error();
	}
      
      }// end switch
  }

} // namespace GRINS
