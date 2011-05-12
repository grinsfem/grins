//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - a low Mach number Navier-Stokes Finite-Element Solver
//
// Copyright (C) 2010,2011 The PECOS Development Team
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
// Definitions for the GRINS::MeshManager class.
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "grins_mesh_manager.h"

#include <iostream>

GRINS::MeshManager::MeshManager()
  : _mesh(NULL),
    _mesh_created_locally(false)
{
  return;
}

GRINS::MeshManager::~MeshManager()
{
  // Only delete the mesh if we actually new'ed it.
  if( this->_mesh_created_locally )
    {
      delete this->_mesh;
      this->_mesh = NULL;
    }
  return;
}

void GRINS::MeshManager::read_input_options( const GetPot& input )
{
  this->_mesh_option = (MESH_OPTION_ENUM)input("mesh-options/mesh_option",
                                               (int)READ_MESH_FROM_FILE);
  this->_print_mesh_info_flag = input("mesh-options/print_mesh_info_flag",
                                      false);

  if(this->_mesh_option==READ_MESH_FROM_FILE)
    {
      if(input.have_variable("mesh-options/mesh_filename"))
        {
          this->_mesh_filename = input("mesh-options/mesh_filename", "NULL");
        }
      else
        {
	  // TODO: Need more consistent error handling.
          std::cerr << " GRINS::MeshManager::read_input_options :" << 
                       " mesh-options/mesh_filename NOT specified " <<
                       std::endl;
          exit(1);
        }
    }
  else if(this->_mesh_option!=MESH_ALREADY_LOADED)
    {
      std::string domain_type_default_value;
      switch (this->_mesh_option)
      {
        case CREATE_1D_MESH:
          domain_type_default_value = "line";
          break;
        case CREATE_2D_MESH:
          domain_type_default_value = "rectangle";
          break;
        case CREATE_3D_MESH:
          domain_type_default_value = "box";
          break;
      }

      this->_domain_type = input("mesh-options/domain_type",
                                 domain_type_default_value);

      this->_domain_x1_min = input("mesh-options/domain_x1_min", 0.0);
      this->_domain_x2_min = input("mesh-options/domain_x2_min", 0.0);
      this->_domain_x3_min = input("mesh-options/domain_x3_min", 0.0);

      this->_domain_x1_max = input("mesh-options/domain_x1_max", 1.0);
      this->_domain_x2_max = input("mesh-options/domain_x2_max", 1.0);
      this->_domain_x3_max = input("mesh-options/domain_x3_max", 1.0);

      this->_mesh_nx1 = input("mesh-options/mesh_nx1", 10);
      this->_mesh_nx2 = input("mesh-options/mesh_nx2", 10);
      this->_mesh_nx3 = input("mesh-options/mesh_nx3", 10);

      // set default element type as libMeshEnums::INVALID_ELEM and
      // switch it to an appropriate one based on domain type
      this->_element_type =
             (libMeshEnums::ElemType)input("mesh-options/element_type",
                                           (int)libMeshEnums::INVALID_ELEM);
    }

  return;
}

libMesh::Mesh* GRINS::MeshManager::get_mesh()
{
  // mesh can be available due to either set_mesh() or build_mesh()
  if( !this->_mesh )
    {
      // TODO: Need more consistent error handling.
      std::cerr << " GRINS::MeshManager::get_mesh :" 
		<< " mesh not yet constructed. " 
		<< std::endl;
      exit(1);
    }

  return this->_mesh;
}

void GRINS::MeshManager::set_mesh( libMesh::Mesh* mesh )
{
  if( this->_mesh )
    {
      // TODO: Need more consistent error handling.
      std::cerr << " GRINS::MeshManager::set_mesh :" << 
                   " mesh is already set " << std::endl;
      exit(1);
    }

  return;
}

void GRINS::MeshManager::build_mesh()
{
  if( this->_mesh_option==MESH_ALREADY_LOADED || this->_mesh )
    {
      // TODO: Need more consistent error handling.
      std::cerr << " GRINS::MeshManager::build_mesh :" << 
                   " mesh is already loaded or set " << std::endl;
      exit(1);
    }

  // Create Mesh object (defaults to dimension 1).
  // According to Roy Stogner, the only read format
  // that won't properly reset the dimension is gmsh.
  this->_mesh = new libMesh::Mesh();

  switch (this->_mesh_option)
  {
    case READ_MESH_FROM_FILE:
      {
        // FIXME: Make this gmsh safe ---  GRINS should worry about this
        (this->_mesh)->read(this->_mesh_filename);
      }
      break;
    case CREATE_1D_MESH:
      {
        if(this->_element_type==libMeshEnums::INVALID_ELEM)
	  {
	    this->_element_type = libMeshEnums::EDGE3;
	  }
	
        libMesh::MeshTools::Generation::build_line(*(this->_mesh),
						   this->_mesh_nx1,
						   this->_domain_x1_min,
						   this->_domain_x1_max,
						   this->_element_type);
	_mesh_created_locally = true;
      }
      break;
    case CREATE_2D_MESH:
      {
        if(this->_element_type==libMeshEnums::INVALID_ELEM)
	  {
	    this->_element_type = libMeshEnums::TRI6;
	  }

	// Reset mesh dimension to 2.
	(this->_mesh)->set_mesh_dimension(2);

        libMesh::MeshTools::Generation::build_square(*(this->_mesh),
						     this->_mesh_nx1,
						     this->_mesh_nx2,
						     this->_domain_x1_min,
						     this->_domain_x1_max,
						     this->_domain_x2_min,
						     this->_domain_x2_max,
						     this->_element_type);
	_mesh_created_locally = true;
      }
      break;
    case CREATE_3D_MESH:
      {
        if(this->_element_type==libMeshEnums::INVALID_ELEM)
	  {
	    this->_element_type = libMeshEnums::TET10;
	  }

	// Reset mesh dimension to 3.
	(this->_mesh)->set_mesh_dimension(3);

        libMesh::MeshTools::Generation::build_cube(*(this->_mesh),
						   this->_mesh_nx1,
						   this->_mesh_nx2,
						   this->_mesh_nx3,
						   this->_domain_x1_min,
						   this->_domain_x1_max,
						   this->_domain_x2_min,
						   this->_domain_x2_max,
						   this->_domain_x3_min,
						   this->_domain_x3_max,
						   this->_element_type);
	_mesh_created_locally = true;
      }
      break;
    default:
      {
	// TODO: Need more consistent error handling.
        std::cerr << " GRINS::MeshManager::build_mesh :" << 
                     " specified mesh-options/mesh_option NOT supported " <<
                     std::endl;
        exit(1);
      }
      break;
  }

  return;
}
