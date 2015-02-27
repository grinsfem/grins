//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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


// C++
#include <iostream>

// This class
#include "grins/grins_enums.h"
#include "grins/mesh_builder.h"

// GRINS
#include "grins/common.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/parsed_function.h"
#include "libmesh/serial_mesh.h"


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

  std::tr1::shared_ptr<libMesh::UnstructuredMesh> MeshBuilder::build
    (const GetPot& input,
     const libMesh::Parallel::Communicator &comm)
  {
    // User needs to tell us if we are generating or reading a mesh
    if( !input.have_variable("mesh-options/mesh_option") &&
        !input.have_variable("MeshOptions/type") )
      {
        libmesh_error_msg("ERROR: Must specify MeshOptions/type in input.");
      }

    // Are we generating the mesh or are we reading one in from a file?
    std::string mesh_build_type = input("MeshOptions/type", "DIE!");
    if( input.have_variable("mesh-options/mesh_option") )
      {
        std::string warning = "WARNING: mesh-options/mesh_option is DEPRECATED.\n";
        warning += "         Please update to use MeshOptions/type.\n";
        warning += "         MeshOptions/type can take values: generate, read\n";
        grins_warning(warning);

        mesh_build_type = input("mesh-options/mesh_option", "DIE!");
      }

    // Make sure the user gave a valid option
    /*! \todo Can remove last 4 checks once mesh-options/mesh_option support is removed. */
    if( mesh_build_type != std::string("generate") &&
        mesh_build_type != std::string("read") &&
        mesh_build_type != std::string("read_mesh_from_file") &&
        mesh_build_type != std::string("create_1D_mesh") &&
        mesh_build_type != std::string("create_2D_mesh") &&
        mesh_build_type != std::string("create_3D_mesh") )
      {
        std::string error = "ERROR: Invalid value of "+mesh_build_type+" for MeshOptions/type.\n";
          error += "       Valid values are: generate\n";
          error += "                         read\n";
          libmesh_error_msg(error);
      }

    // Create UnstructuredMesh object (defaults to dimension 1).
    libMesh::UnstructuredMesh* mesh;

    // Were we specifically asked to use a ParallelMesh or SerialMesh?
    {
      std::string mesh_class = input("MeshOptions/mesh_class", "default");

      if( input.have_variable("mesh-options/mesh_class") )
        {
          std::string warning = "WARNING: mesh-options/mesh_class is DEPRECATED.\n";
          warning += "         Please update to use MeshOptions/mesh_class.\n";
          warning += "         MeshOptions/mesh_class can take values: serial, parallel\n";
          grins_warning(warning);

          mesh_class = input("mesh-options/mesh_class", "default");
        }

      if (mesh_class == "parallel")
        mesh = new libMesh::ParallelMesh(comm);
      else if (mesh_class == "serial")
        mesh = new libMesh::SerialMesh(comm);
      else if (mesh_class == "default")
        mesh = new libMesh::Mesh(comm);
      else
        {
          std::string error = "ERROR: Invalid mesh_class "+mesh_class+" input for MeshOptions/mesh_class.\n";
            error += "       Valid choices are: serial, parallel.\n";
            libmesh_error_msg(error);
        }
    }

    // Read mesh from file
    if(mesh_build_type =="read_mesh_from_file" /* This is deprecated */ ||
       mesh_build_type == "read" )
      {
        // Make sure the user set the filename to read
        if( !input.have_variable("mesh-options/mesh_filename") /* This is deprecated */ &&
            !input.have_variable("MeshOptions/Read/filename") )
          {
            libmesh_error_msg("ERROR: Must specify MeshOptions/Read/filename for reading mesh.");
          }

        std::string mesh_filename = input("MeshOptions/Read/filename", "DIE!");

        if( input.have_variable("mesh-options/mesh_filename") )
          {
            std::string warning = "WARNING: mesh-options/mesh_filename is DEPRECATED.\n";
            warning += "         Please update to use MeshOptions/Read/filename.\n";
            grins_warning(warning);

            mesh_filename = input("mesh-options/mesh_filename", "DIE!");
          }

	// According to Roy Stogner, the only read format
	// that won't properly reset the dimension is gmsh.
	/*! \todo Need to a check a GMSH meshes */
	mesh->read(mesh_filename);
      }

    // Generate the mesh using built-in libMesh functions
    else if(mesh_build_type=="create_1D_mesh" /* This is deprecated */ ||
            mesh_build_type=="create_2D_mesh" /* This is deprecated */ ||
            mesh_build_type=="create_3D_mesh" /* This is deprecated */ ||
            mesh_build_type=="generate")
      {
        this->generate_mesh(mesh_build_type,input,mesh);
      }

    else
      {
        // Shouldn't have gotten here
        libmesh_error();
      }

    /* Only do the mesh refinement here if we don't have a restart file.
       Otherwise, we need to wait until we've read in the restart file.
       That is done in Simulation::check_for_restart */
    if( !input.have_variable("restart-options/restart_file") )
      {
        this->do_mesh_refinement_from_input( input, comm, *mesh );
      }

    return std::tr1::shared_ptr<libMesh::UnstructuredMesh>(mesh);
  }

  void MeshBuilder::generate_mesh( const std::string& mesh_build_type, const GetPot& input,
                                   libMesh::UnstructuredMesh* mesh )
  {
    unsigned int dimension = input("MeshOptions/Generation/dimension",0);

    if( !input.have_variable("mesh-options/mesh_option") /* Deprecated */ &&
        !input.have_variable("MeshOptions/Generation/dimension") )
      {
        libmesh_error_msg("ERROR: Must specify MeshOptions/Generation/dimension for generating mesh.");
      }

    /* Remove these once suport for mesh-options/mesh_option is removed */
    if( mesh_build_type == "create_1D_mesh" /* This is deprecated */ )
      dimension = 1;

    if( mesh_build_type=="create_2D_mesh" /* This is deprecated */ )
      dimension = 2;

    if( mesh_build_type=="create_3D_mesh" /* This is deprecated */ )
      dimension = 3;

    // Set the mesh dimension
    mesh->set_mesh_dimension(dimension);

    libMesh::Real x_min = input("MeshOptions/Generation/x_min", 0.0);
    if( input.have_variable("mesh-options/domain_x1_min") )
      {
        std::string warning = "WARNING: mesh-options/domain_x1_min is DEPRECATED.\n";
          warning += "         Please update to use MeshOptions/Generation/x_min.\n";
          grins_warning(warning);

          x_min = input("mesh-options/domain_x1_min", 0.0);
      }

    libMesh::Real y_min = input("MeshOptions/Generation/y_min", 0.0);
    if( input.have_variable("mesh-options/domain_x2_min") )
      {
        std::string warning = "WARNING: mesh-options/domain_x2_min is DEPRECATED.\n";
          warning += "         Please update to use MeshOptions/Generation/y_min.\n";
          grins_warning(warning);

          y_min = input("mesh-options/domain_x2_min", 0.0);
      }

    libMesh::Real z_min = input("MeshOptions/Generation/z_min", 0.0);
    if( input.have_variable("mesh-options/domain_x3_min") )
      {
        std::string warning = "WARNING: mesh-options/domain_x3_min is DEPRECATED.\n";
          warning += "         Please update to use MeshOptions/Generation/z_min.\n";
          grins_warning(warning);

          z_min = input("mesh-options/domain_x3_min", 0.0);
      }

    libMesh::Real x_max = input("MeshOptions/Generation/x_max", 1.0);
    if( input.have_variable("mesh-options/domain_x1_max") )
      {
        std::string warning = "WARNING: mesh-options/domain_x1_max is DEPRECATED.\n";
          warning += "         Please update to use MeshOptions/Generation/x_max.\n";
          grins_warning(warning);

          x_max = input("mesh-options/domain_x1_max", 1.0);
      }

    libMesh::Real y_max = input("MeshOptions/Generation/y_max", 1.0);
    if( input.have_variable("mesh-options/domain_x2_max") )
      {
        std::string warning = "WARNING: mesh-options/domain_x2_max is DEPRECATED.\n";
          warning += "         Please update to use MeshOptions/Generation/y_max.\n";
          grins_warning(warning);

          y_max = input("mesh-options/domain_x2_max", 1.0);
      }

    libMesh::Real z_max = input("MeshOptions/Generation/z_max", 1.0);
    if( input.have_variable("mesh-options/domain_x3_max") )
      {
        std::string warning = "WARNING: mesh-options/domain_x3_max is DEPRECATED.\n";
          warning += "         Please update to use MeshOptions/Generation/z_max.\n";
          grins_warning(warning);

          z_max = input("mesh-options/domain_x3_max", 1.0);
      }

    // Make sure user gave us info about how many elements to use
    if( !input.have_variable("mesh-options/mesh_nx1") /* Deprecated */ &&
        !input.have_variable("MeshOptions/Generation/n_elems_x") )
      {
        libmesh_error_msg("ERROR: Must supply MeshOptions/Generation/n_elems_x for mesh generation.");
      }

    if( dimension > 1 )
      {
        if( !input.have_variable("mesh-options/mesh_nx2") /* Deprecated */ &&
        !input.have_variable("MeshOptions/Generation/n_elems_y") )
          {
            libmesh_error_msg("ERROR: Must supply MeshOptions/Generation/n_elems_y for mesh generation.");
          }
      }

    if( dimension > 2 )
      {
        if( !input.have_variable("mesh-options/mesh_nx3") /* Deprecated */ &&
        !input.have_variable("MeshOptions/Generation/n_elems_z") )
          {
            libmesh_error_msg("ERROR: Must supply MeshOptions/Generation/n_elems_z for mesh generation.");
          }
      }

    unsigned int n_elems_x = input("MeshOptions/Generation/n_elems_x", 0);
    if( input.have_variable("mesh-options/mesh_nx1") )
      {
        std::string warning = "WARNING: mesh-options/mesh_nx1 is DEPRECATED.\n";
        warning += "         Please update to use MeshOptions/Generation/n_elems_x.\n";
        grins_warning(warning);

        n_elems_x = input("mesh-options/mesh_nx1", 0);
      }

    unsigned int n_elems_y = input("MeshOptions/Generation/n_elems_y", 0);
    if( input.have_variable("mesh-options/mesh_nx2") )
      {
        std::string warning = "WARNING: mesh-options/mesh_nx2 is DEPRECATED.\n";
        warning += "         Please update to use MeshOptions/Generation/n_elems_y.\n";
        grins_warning(warning);

        n_elems_y = input("mesh-options/mesh_nx2", 0);
      }

    unsigned int n_elems_z = input("MeshOptions/Generation/n_elems_z", 0);
    if( input.have_variable("mesh-options/mesh_nx3") )
      {
        std::string warning = "WARNING: mesh-options/mesh_nx3 is DEPRECATED.\n";
        warning += "         Please update to use MeshOptions/Generation/n_elems_z.\n";
        grins_warning(warning);

        n_elems_z = input("mesh-options/mesh_nx3", 0);
      }

    std::string element_type = input("MeshOptions/Generation/element_type", "default");
    if( input.have_variable("mesh-options/element_type") )
      {
        std::string warning = "WARNING: mesh-options/element_type is DEPRECATED.\n";
        warning += "         Please update to use MeshOptions/Generation/element_type.\n";
        grins_warning(warning);

        element_type = input("mesh-options/element_type", "default");
      }

    if( dimension == 1 )
      {
        if(element_type=="default")
	  {
	    element_type = "EDGE3";
	  }

        GRINSEnums::ElemType element_enum_type =
	  libMesh::Utility::string_to_enum<GRINSEnums::ElemType>(element_type);

	libMesh::MeshTools::Generation::build_line(*mesh,
                                                   n_elems_x,
                                                   x_min,
                                                   x_max,
						   element_enum_type);
      }

    else if( dimension == 2 )
      {
	if(element_type=="default")
	  {
	    element_type = "TRI6";
	  }

        GRINSEnums::ElemType element_enum_type =
	  libMesh::Utility::string_to_enum<GRINSEnums::ElemType>(element_type);

	libMesh::MeshTools::Generation::build_square(*mesh,
                                                     n_elems_x,
                                                     n_elems_y,
                                                     x_min,
                                                     x_max,
                                                     y_min,
                                                     y_max,
						     element_enum_type);
      }

    else if( dimension == 3 )
      {
	if(element_type=="default")
	  {
	    element_type = "TET10";
	  }

        GRINSEnums::ElemType element_enum_type =
	  libMesh::Utility::string_to_enum<GRINSEnums::ElemType>(element_type);

	libMesh::MeshTools::Generation::build_cube(*mesh,
                                                   n_elems_x,
                                                   n_elems_y,
                                                   n_elems_z,
                                                   x_min,
                                                   x_max,
                                                   y_min,
                                                   y_max,
                                                   z_min,
                                                   z_max,
						   element_enum_type);
      }

    else
      {
        // This shouldn't have happened
	libmesh_error();
      }
  }

  void MeshBuilder::do_mesh_refinement_from_input( const GetPot& input,
                                                   const libMesh::Parallel::Communicator &comm,
                                                   libMesh::UnstructuredMesh& mesh ) const
  {
    std::string redistribution_function_string =
            input("mesh-options/redistribute", std::string("0"));

    if (redistribution_function_string != "0")
      {
        libMesh::ParsedFunction<libMesh::Real>
          redistribution_function(redistribution_function_string);

        libMesh::MeshTools::Modification::redistribute
          (mesh, redistribution_function);

        // Redistribution can create distortions *within* second-order
        // elements, which can then be magnified by refinement.  Let's
        // undistort everything by converting to first order and back
        // if necessary.

        // FIXME - this only works for meshes with uniform geometry
        // order equal to FIRST or (full-order) SECOND.

        const libMesh::Elem *elem = *mesh.elements_begin();

        if (elem->default_order() != libMesh::FIRST)
          {
            mesh.all_first_order();
            mesh.all_second_order();
          }
      }

    int uniformly_refine = input("mesh-options/uniformly_refine", 0);
    
    if( uniformly_refine > 0 )
      {
        libMesh::MeshRefinement(mesh).uniformly_refine(uniformly_refine);
      }

    std::string h_refinement_function_string =
            input("mesh-options/locally_h_refine", std::string("0"));

    if (h_refinement_function_string != "0")
      { 
        libMesh::ParsedFunction<libMesh::Real>
          h_refinement_function(h_refinement_function_string);

        libMesh::MeshRefinement mesh_refinement(mesh);

        libMesh::dof_id_type found_refinements = 0;
        do {
          found_refinements = 0;
          unsigned int max_level_refining = 0;

          libMesh::MeshBase::element_iterator elem_it =
                  mesh.active_elements_begin();
          libMesh::MeshBase::element_iterator elem_end =
                  mesh.active_elements_end();
          for (; elem_it != elem_end; ++elem_it)
            {
              libMesh::Elem *elem = *elem_it;

              const libMesh::Real refinement_val =
                h_refinement_function(elem->centroid());

	      const unsigned int n_refinements = refinement_val > 0 ?
                refinement_val : 0;

              if (elem->level() - uniformly_refine < n_refinements)
                {
                  elem->set_refinement_flag(libMesh::Elem::REFINE);
                  found_refinements++;
		  max_level_refining = std::max(max_level_refining,
                                                elem->level());
                }
            }

          comm.max(found_refinements);
          comm.max(max_level_refining);

          if (found_refinements)
            {
	      std::cout << "Found up to " << found_refinements << 
                " elements to refine on each processor," << std::endl;
	      std::cout << "with max level " << max_level_refining << std::endl;
              mesh_refinement.refine_and_coarsen_elements();

              if( input.have_variable("restart-options/restart_file") )
                {
                  std::cout << "Warning: it is known that locally_h_refine is broken when restarting." << std::endl
                            << "         and multiple refinement passes are done. We are forcibly" << std::endl
                            << "         limiting the refinement to one pass until this issue is resolved." << std::endl;
                  break;
                }
            }

        } while(found_refinements);

      }

    return;
  }

} // namespace GRINS
