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

namespace GRINS
{

  MagnetostaticsBCHandling::MagnetostaticsBCHandling( const std::string& physics_name,
						      const GetPot& input)
    : BCHandlingBase(physics_name)
  {
    _A_var_name = input("Physics/VariableNames/MagneticPotential", GRINS::A_var_name_default );
    
    std::string id_str = "Physics/"+_physics_name+"/bc_ids";
    std::string bc_str = "Physics/"+_physics_name+"/bc_types";
    
    this->read_bc_data( input, id_str, bc_str );
    
    return;
  }

  MagnetostaticsBCHandling::~MagnetostaticsBCHandling()
  {
    return;
  }

  int MagnetostaticsBCHandling::string_to_int( const std::string& bc_type ) const
  {
    AM_BC_TYPES bc_type_out;
    
    if( bc_type == "axisymmetric" )
      bc_type_out = AXISYMMETRIC;
    
    else if( bc_type == "unbounded" )
      bc_type_out = UNBOUNDED;

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
  
  void MagnetostaticsBCHandling::init_bc_data( const BoundaryID bc_id, 
					       const std::string& bc_id_string, 
					       const int bc_type, 
					       const GetPot& input )
  {
    switch(bc_type)
      {
      case(AXISYMMETRIC):
	{
	
	}
	break;
      case(UNBOUNDED):
	// Do nothing BC
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

  void MagnetostaticsBCHandling::user_apply_neumann_bcs( libMesh::FEMContext& context,
							 VariableIndex var,
							 bool request_jacobian,
							 BoundaryID bc_id,
							 BCType bc_type ) const
  {
    switch( bc_type )
      {
      
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
