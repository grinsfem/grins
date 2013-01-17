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

#ifndef MESH_BUILDER_H
#define MESH_BUILDER_H

#include "boost/tr1/memory.hpp"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/libmesh.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_refinement.h"

namespace GRINS
{

  class MeshBuilder
  {    
  public:

    //! This Object handles building a libMesh::Mesh
    /*! Based on runtime input, either a generic 1, 2, or 3-dimensional
        mesh is built; or is read from input from a specified file. */
    MeshBuilder();
    ~MeshBuilder();

    void read_input_options( const GetPot& input );

    //! Builds the libMesh::Mesh according to input options.
    std::tr1::shared_ptr<libMesh::Mesh> build(const GetPot& input );

  };

} //End namespace block

#endif //MESH_BUILDER_H
