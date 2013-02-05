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

// C++
#include <iostream>

// This class
#include "grins/mesh_builder.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_refinement.h"


namespace GRINS
{

  MeshBuilder::MeshBuilder()
  {
    return;
  }

  MeshBuilder::~MeshBuilder()
  {
    return;
  }

  std::tr1::shared_ptr<libMesh::Mesh> MeshBuilder::build(const GetPot& input)
  {
    // First, read all needed variables
    std::string mesh_option = input("mesh-options/mesh_option", "NULL");
    std::string mesh_filename = input("mesh-options/mesh_filename", "NULL");

    libMesh::Real domain_x1_min = input("mesh-options/domain_x1_min", 0.0);
    libMesh::Real domain_x2_min = input("mesh-options/domain_x2_min", 0.0);
    libMesh::Real domain_x3_min = input("mesh-options/domain_x3_min", 0.0);

    libMesh::Real domain_x1_max = input("mesh-options/domain_x1_max", 1.0); 
    libMesh::Real domain_x2_max = input("mesh-options/domain_x2_max", 1.0);
    libMesh::Real domain_x3_max = input("mesh-options/domain_x3_max", 1.0);

    int mesh_nx1 = input("mesh-options/mesh_nx1", -1);
    int mesh_nx2 = input("mesh-options/mesh_nx2", -1);
    int mesh_nx3 = input("mesh-options/mesh_nx3", -1);

    int uniformly_refine = input("mesh-options/uniformly_refine", 0);
    
    std::string element_type = input("mesh-options/element_type", "NULL");

    // Make sure the user told us what to do
    if(mesh_option == "NULL")
      {
	std::cerr << " MeshBuilder::read_input_options :"
		  << " mesh-options/mesh_option NOT specified "
		  << std::endl;
	libmesh_error();
      }

    // Create Mesh object (defaults to dimension 1).
    libMesh::Mesh* mesh = new libMesh::Mesh();

    if(mesh_option=="read_mesh_from_file")
      {
	// According to Roy Stogner, the only read format
	// that won't properly reset the dimension is gmsh.
	/*! \todo Need to a check a GMSH meshes */
	mesh->read(mesh_filename);
      }

    else if(mesh_option=="create_1D_mesh")
      {
	if(element_type=="NULL")
	  {
	    element_type = "EDGE3";
	  }
      
	libMeshEnums::ElemType element_enum_type =
	  libMesh::Utility::string_to_enum<libMeshEnums::ElemType>(element_type);
      
	libMesh::MeshTools::Generation::build_line(*mesh,
						   mesh_nx1,
						   domain_x1_min,
						   domain_x1_max,
						   element_enum_type);
      }
      
    else if(mesh_option=="create_2D_mesh")
      {
	if(element_type=="NULL")
	  {
	    element_type = "TRI6";
	  }

	// Reset mesh dimension to 2.
	mesh->set_mesh_dimension(2);

	libMeshEnums::ElemType element_enum_type =
	  libMesh::Utility::string_to_enum<libMeshEnums::ElemType>(element_type);

	libMesh::MeshTools::Generation::build_square(*mesh,
						     mesh_nx1,
						     mesh_nx2,
						     domain_x1_min,
						     domain_x1_max,
						     domain_x2_min,
						     domain_x2_max,
						     element_enum_type);
      }

    else if(mesh_option=="create_3D_mesh")
      {
	if(element_type=="NULL")
	  {
	    element_type = "TET10";
	  }

	// Reset mesh dimension to 3.
	mesh->set_mesh_dimension(3);

	libMeshEnums::ElemType element_enum_type =
	  libMesh::Utility::string_to_enum<libMeshEnums::ElemType>(element_type);

	libMesh::MeshTools::Generation::build_cube(*mesh,
						   mesh_nx1,
						   mesh_nx2,
						   mesh_nx3,
						   domain_x1_min,
						   domain_x1_max,
						   domain_x2_min,
						   domain_x2_max,
						   domain_x3_min,
						   domain_x3_max,
						   element_enum_type);
      }

    else
      {
	std::cerr << " MeshBuilder::build_mesh :"
		  << " mesh-options/mesh_option [" << mesh_option
		  << "] NOT supported " << std::endl;
	libmesh_error();
      }

    if( uniformly_refine > 0 )
      {
	libMesh::MeshRefinement(*mesh).uniformly_refine(uniformly_refine);
      }

    return std::tr1::shared_ptr<libMesh::Mesh>(mesh);
  }

} // namespace GRINS
