//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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

// libMesh
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/partitioner.h"
#include "libmesh/parsed_function.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/elem.h"
#include "libmesh/enum_order.h"

namespace GRINS
{
  std::shared_ptr<libMesh::UnstructuredMesh> MeshBuilder::build
  ( const GetPot & input,
    const libMesh::Parallel::Communicator & comm )
  {
    // User needs to tell us if we are generating or reading a mesh
    // We infer this by checking and seeing if the use has a Mesh/Read
    // or a Mesh/Generation section
    if( !input.have_section("Mesh/Read/") &&
        !input.have_section("Mesh/Generation/") )
      libmesh_error_msg("ERROR: Must specify either Mesh/Read or Mesh/Generation in input.");

    // But you can't have it both ways
    if( input.have_section("Mesh/Read/") &&
        input.have_section("Mesh/Generation/") )
      libmesh_error_msg("ERROR: Can only specify one of Mesh/Read and Mesh/Generation");

    // Are we generating the mesh or are we reading one in from a file?
    std::string mesh_build_type = "NULL";
    if( input.have_section("Mesh/Read/") )
      mesh_build_type = "read";

    else if( input.have_section("Mesh/Generation/") )
      mesh_build_type = "generate";

    // Make sure the user gave a valid option
    /*! \todo Can remove last 4 checks once mesh-options/mesh_option support is removed. */
    if( mesh_build_type != std::string("generate") &&
        mesh_build_type != std::string("read") )
      {
        std::string error = "ERROR: Invalid value of "+mesh_build_type+" for Mesh/type.\n";
        error += "       Valid values are: generate\n";
        error += "                         read\n";
        libmesh_error_msg(error);
      }

    // Create UnstructuredMesh object (defaults to dimension 1).
    std::shared_ptr<libMesh::UnstructuredMesh> mesh;

    // Were we specifically asked to use a ParallelMesh or SerialMesh?
    {
      std::string mesh_class = input("Mesh/class", "default");

      if (mesh_class == "parallel")
        mesh = std::make_shared<libMesh::ParallelMesh>(comm);
      else if (mesh_class == "serial")
        mesh = std::make_shared<libMesh::SerialMesh>(comm);
      else if (mesh_class == "default")
        mesh = std::make_shared<libMesh::Mesh>(comm);
      else
        {
          std::string error = "ERROR: Invalid class "+mesh_class+" input for Mesh/class.\n";
          error += "       Valid choices are: serial, parallel.\n";
          libmesh_error_msg(error);
        }
    }

    // Read mesh from file
    if( mesh_build_type == "read" )
      {
        // Make sure the user set the filename to read
        if( !input.have_variable("Mesh/Read/filename") )
          libmesh_error_msg("ERROR: Must specify Mesh/Read/filename for reading mesh.");

        std::string mesh_filename = input("Mesh/Read/filename", "DIE!");

        // According to Roy Stogner, the only read format
        // that won't properly reset the dimension is gmsh.
        /*! \todo Need to a check a GMSH meshes */
        mesh->read(mesh_filename);

        // If we have a first order mesh file but we need second order
        // elements we should fix that.
        bool all_second_order = input("Mesh/all_second_order", false);
        if (all_second_order)
          mesh->all_second_order();
      }

    // Generate the mesh using built-in libMesh functions
    else if(mesh_build_type=="generate")
      this->generate_mesh(input,*mesh);

    // Shouldn't have gotten here
    else
      libmesh_error();

    /* Only do the mesh refinement here if we don't have a restart file.
       Otherwise, we need to wait until we've read in the restart file.
       That is done in Simulation::check_for_restart */
    if( !input.have_variable("restart-options/restart_file") )
      this->do_mesh_refinement_from_input( input, comm, *mesh );

    return mesh;
  }

  void MeshBuilder::generate_mesh( const GetPot & input,
                                   libMesh::UnstructuredMesh & mesh )
  {
    unsigned int dimension = input("Mesh/Generation/dimension",0);

    if( !input.have_variable("Mesh/Generation/dimension") )
      libmesh_error_msg("ERROR: Must specify Mesh/Generation/dimension for generating mesh.");

    // Set the mesh dimension
    mesh.set_mesh_dimension(dimension);

    /* Now look for spatial extent of the grid that the user wants to generate. */
    libMesh::Real x_min = input("Mesh/Generation/x_min", 0.0);
    libMesh::Real x_max = input("Mesh/Generation/x_max", 1.0);

    /* We only check the y_{min,max} input if dimension is > 1 so that GetPot
       UFO detection will give us an error if we have this in the input file
       and are only using a 1D grid. */
    libMesh::Real y_min = 0.0;
    libMesh::Real y_max = 0.0;

    if( dimension > 1 )
      {
        y_min = input("Mesh/Generation/y_min", 0.0);
        y_max = input("Mesh/Generation/y_max", 1.0);
      }

    /* We only check the z_{min,max} input if dimension is > 2 so that GetPot
       UFO detection will give us an error if we have this in the input file
       and are only using a 1D or 2D grid. */
    libMesh::Real z_min = 0.0;
    libMesh::Real z_max = 0.0;

    if( dimension > 2 )
      {
        z_min = input("Mesh/Generation/z_min", 0.0);
        z_max = input("Mesh/Generation/z_max", 1.0);
      }

    /* Now check for the number of elements in each direction */

    // Make sure user gave us info about how many elements to use
    if( !input.have_variable("Mesh/Generation/n_elems_x") )
      libmesh_error_msg("ERROR: Must supply Mesh/Generation/n_elems_x for mesh generation.");

    unsigned int n_elems_x = input("Mesh/Generation/n_elems_x", 0);

    /* We only check n_elems_y input if dimension is > 1 so that GetPot
       UFO detection will give us an error if we have this in the input file
       and are only using a 1D grid. */
    unsigned int n_elems_y = 0;
    if( dimension > 1 )
      {
        if( !input.have_variable("Mesh/Generation/n_elems_y") )
          libmesh_error_msg("ERROR: Must supply Mesh/Generation/n_elems_y for mesh generation.");

        n_elems_y = input("Mesh/Generation/n_elems_y", 0);
      }

    /* We only check n_elems_z input if dimension is > 2 so that GetPot
       UFO detection will give us an error if we have this in the input file
       and are only using a 1D or 2D grid. */
    unsigned int n_elems_z = 0;
    if( dimension > 2 )
      {
        if( !input.have_variable("Mesh/Generation/n_elems_z") )
          libmesh_error_msg("ERROR: Must supply Mesh/Generation/n_elems_z for mesh generation.");

        n_elems_z = input("Mesh/Generation/n_elems_z", 0);
      }

    /* Now grab the element_type the user wants for the mesh. */

    std::string element_type = input("Mesh/Generation/element_type", "default");

    /* Now generate the mesh. */
    if( dimension == 1 )
      {
        if(element_type=="default")
          element_type = "EDGE3";

        GRINSEnums::ElemType element_enum_type =
          libMesh::Utility::string_to_enum<GRINSEnums::ElemType>(element_type);

        libMesh::MeshTools::Generation::build_line(mesh,
                                                   n_elems_x,
                                                   x_min,
                                                   x_max,
                                                   element_enum_type);
      }

    else if( dimension == 2 )
      {
        if(element_type=="default")
          element_type = "TRI6";

        GRINSEnums::ElemType element_enum_type =
          libMesh::Utility::string_to_enum<GRINSEnums::ElemType>(element_type);

        libMesh::MeshTools::Generation::build_square(mesh,
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
          element_type = "TET10";

        GRINSEnums::ElemType element_enum_type =
          libMesh::Utility::string_to_enum<GRINSEnums::ElemType>(element_type);

        libMesh::MeshTools::Generation::build_cube(mesh,
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

  void MeshBuilder::do_mesh_refinement_from_input( const GetPot & input,
                                                   const libMesh::Parallel::Communicator & comm,
                                                   libMesh::UnstructuredMesh & mesh ) const
  {
    std::string redistribution_function_string =
      input("Mesh/Redistribution/function", std::string("0"));

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

    bool allow_remote_elem_deletion = input("Mesh/Refinement/allow_remote_elem_deletion", true);
    if (allow_remote_elem_deletion == false)
      mesh.allow_remote_element_removal(false);

    bool allow_renumbering = input("Mesh/Refinement/allow_renumbering", true);
    if (allow_renumbering == false)
      mesh.allow_renumbering(false);

    bool disable_partitioning = input("Mesh/Refinement/disable_partitioning", false);
    if (disable_partitioning == true)
      mesh.partitioner() = nullptr;

    int uniformly_refine = input("Mesh/Refinement/uniformly_refine", 0);

    if( uniformly_refine > 0 )
      libMesh::MeshRefinement(mesh).uniformly_refine(uniformly_refine);

    std::string h_refinement_function_string =
      input("Mesh/Refinement/locally_h_refine", std::string("0"));

    if (h_refinement_function_string != "0")
      {
        libMesh::ParsedFunction<libMesh::Real>
          h_refinement_function(h_refinement_function_string);

        libMesh::MeshRefinement mesh_refinement(mesh);

        libMesh::dof_id_type found_refinements = 0;

        do {
          found_refinements = 0;
          unsigned int max_level_refining = 0;

          for (auto & elem : mesh.active_element_ptr_range())
            {
              const libMesh::Real refinement_val =
                h_refinement_function(elem->vertex_average());

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

      } // if (h_refinement_function_string != "0")
  }

} // namespace GRINS
