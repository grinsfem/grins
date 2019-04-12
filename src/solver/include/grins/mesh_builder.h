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


#ifndef GRINS_MESH_BUILDER_H
#define GRINS_MESH_BUILDER_H

// libMesh
#include "libmesh/mesh.h"
#include "libmesh/getpot.h"


// GRINS
#include "grins/common.h"
#include <memory>

namespace GRINS
{

  class MeshBuilder
  {
  public:

    //! This Object handles building a libMesh::UnstructuredMesh subclass.
    /*! Based on runtime input, either a generic 1, 2, or 3-dimensional
      mesh is built; or is read from input from a specified file. */
    MeshBuilder(){};
    ~MeshBuilder(){};

    //! Builds the libMesh::Mesh according to input options.
    std::shared_ptr<libMesh::UnstructuredMesh> build
    ( const GetPot& input,
      const libMesh::Parallel::Communicator &comm
      LIBMESH_CAN_DEFAULT_TO_COMMWORLD );

    //! Refine the mesh based on user input parameters
    /*! There are several parameters that allow for the user to specify
      mesh refinements. They are factored out here because where we
      need to apply the refinements depends on the existence of a
      restart file. */
    void do_mesh_refinement_from_input( const GetPot& input,
                                        const libMesh::Parallel::Communicator &comm,
                                        libMesh::UnstructuredMesh& mesh ) const;

  private:

    void generate_mesh( const std::string& mesh_build_type, const GetPot& input,
                        libMesh::UnstructuredMesh* mesh );

    //! Helper function for displaying deprecated warnings.
    template <typename T>
    void deprecated_option( const GetPot& input, const std::string& old_option,
                            const std::string& new_option, const T& default_value,
                            T& option_value ) const;

  };

  template <typename T>
  inline
  void MeshBuilder::deprecated_option( const GetPot& input, const std::string& old_option,
                                       const std::string& new_option, const T& default_value,
                                       T& option_value ) const
  {
    if( input.have_variable(old_option) )
      {
        std::string warning = "WARNING: "+old_option+" is DEPRECATED.\n";
        warning += "         Please update to use "+new_option+".\n";
        grins_warning(warning);

        option_value = input(old_option, default_value);
      }

    return;
  }
} // end namespace block

#endif // GRINS_MESH_BUILDER_H
