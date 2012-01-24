//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
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
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "mesh_builder.h"

#include <iostream>

GRINS::MeshBuilder::MeshBuilder( const GetPot& input )
{
  this->read_input_options( input );
  return;
}

GRINS::MeshBuilder::~MeshBuilder()
{
  return;
}

void GRINS::MeshBuilder::read_input_options( const GetPot& input )
{
  // First check the user told us what to do to generate a mesh
  if(input.have_variable("mesh-options/mesh_option"))
    {
      this->_mesh_option = input("mesh-options/mesh_option", "NULL");
    }
  else
    {
      std::cerr << " GRINS::MeshBuilder::read_input_options :"
                << " mesh-options/mesh_option NOT specified "
                << std::endl;
      libmesh_error();
    }

  // If they want it read from a file, stash the filename
  if(this->_mesh_option=="read_mesh_from_file")
    {
      if(input.have_variable("mesh-options/mesh_filename"))
        {
          this->_mesh_filename = input("mesh-options/mesh_filename", "NULL");
        }
      else
        {
	  // TODO: Need more consistent error handling.
          std::cerr << " GRINS::MeshBuilder::read_input_options :"
                    << " mesh-options/mesh_filename NOT specified "
                    << std::endl;
          libmesh_error();
        }
    }
  // Otherwise, stash the necessary data to create a mesh
  else
    {
      this->_domain_x1_min = input("mesh-options/domain_x1_min", 0.0);
      this->_domain_x2_min = input("mesh-options/domain_x2_min", 0.0);
      this->_domain_x3_min = input("mesh-options/domain_x3_min", 0.0);

      this->_domain_x1_max = input("mesh-options/domain_x1_max", 1.0);
      this->_domain_x2_max = input("mesh-options/domain_x2_max", 1.0);
      this->_domain_x3_max = input("mesh-options/domain_x3_max", 1.0);

      this->_mesh_nx1 = input("mesh-options/mesh_nx1", 10);
      this->_mesh_nx2 = input("mesh-options/mesh_nx2", 10);
      this->_mesh_nx3 = input("mesh-options/mesh_nx3", 10);

      this->_element_type = input("mesh-options/element_type", "NULL");
    }

  return;
}

libMesh::AutoPtr<libMesh::Mesh> GRINS::MeshBuilder::build()
{
  // Create Mesh object (defaults to dimension 1).
  libMesh::Mesh* mesh = new libMesh::Mesh();

  if(this->_mesh_option=="read_mesh_from_file")
    {
      // According to Roy Stogner, the only read format
      // that won't properly reset the dimension is gmsh.
      /*! \todo Need to a check a GMSH meshes */
      mesh->read(this->_mesh_filename);
    }
  else if(this->_mesh_option=="create_1D_mesh")
    {
      if(this->_element_type=="NULL")
	{
	  this->_element_type = "EDGE3";
	}

      libMeshEnums::ElemType _element_enum_type =
                      libMesh::Utility::string_to_enum<libMeshEnums::ElemType>(this->_element_type);

      libMesh::MeshTools::Generation::build_line(*mesh,
						 this->_mesh_nx1,
						 this->_domain_x1_min,
						 this->_domain_x1_max,
						 _element_enum_type);
    }
  else if(this->_mesh_option=="create_2D_mesh")
    {
      if(this->_element_type=="NULL")
	{
	  this->_element_type = "TRI6";
	}

      // Reset mesh dimension to 2.
      mesh->set_mesh_dimension(2);

      libMeshEnums::ElemType _element_enum_type =
                      libMesh::Utility::string_to_enum<libMeshEnums::ElemType>(this->_element_type);

      libMesh::MeshTools::Generation::build_square(*mesh,
						   this->_mesh_nx1,
						   this->_mesh_nx2,
						   this->_domain_x1_min,
						   this->_domain_x1_max,
						   this->_domain_x2_min,
						   this->_domain_x2_max,
						   _element_enum_type);
    }
  else if(this->_mesh_option=="create_3D_mesh")
    {
      if(this->_element_type=="NULL")
	{
	  this->_element_type = "TET10";
	}

      // Reset mesh dimension to 3.
      mesh->set_mesh_dimension(3);

      libMeshEnums::ElemType _element_enum_type =
                      libMesh::Utility::string_to_enum<libMeshEnums::ElemType>(this->_element_type);

      libMesh::MeshTools::Generation::build_cube(*mesh,
						 this->_mesh_nx1,
						 this->_mesh_nx2,
						 this->_mesh_nx3,
						 this->_domain_x1_min,
						 this->_domain_x1_max,
						 this->_domain_x2_min,
						 this->_domain_x2_max,
						 this->_domain_x3_min,
						 this->_domain_x3_max,
						 _element_enum_type);
    }
  else
    {
      std::cerr << " GRINS::MeshBuilder::build_mesh :"
                << " mesh-options/mesh_option [" << this->_mesh_option
                << "] NOT supported " << std::endl;
      libmesh_error();
    }

  return AutoPtr(mesh);
}
