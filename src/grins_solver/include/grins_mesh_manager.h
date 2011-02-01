//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2010 The PECOS Development Team
//
// Please see http://pecos.ices.utexas.edu for more information.
//
// This file is part of GRINS.
//
// GRINS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GRINS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GRINS.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------
//
// Declarations for the GRINS::MeshManager class.
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_MESH_MANAGER_H
#define GRINS_MESH_MANAGER_H

#include "libmesh.h"
#include "mesh.h"
#include "getpot.h"

namespace GRINS
{

  class MeshManager
  {
    
  public:
    // enumeration for mesh options
    enum MESH_OPTION_ENUM { MESH_ALREADY_LOADED = 0,
                            READ_MESH_FROM_FILE = 1,
                            CREATE_1D_MESH      = 2,
                            CREATE_2D_MESH      = 3,
                            CREATE_3D_MESH      = 4 };

    MeshManager();
    ~MeshManager();

    void read_input_options( const GetPot& input );

    // get/set pair for mesh object
    libMesh::Mesh* get_mesh(); 
    void set_mesh( libMesh::Mesh* mesh );

  private:
    int _mesh_option;
    bool _print_mesh_info_flag;

    libMesh::Mesh* _mesh;
  };

} //End namespace block

#endif //GRINS_MESH_MANAGER_H
