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
// Declarations for the GRINS::MeshManager class.
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_MESH_MANAGER_H
#define GRINS_MESH_MANAGER_H

#include "getpot.h"
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"

namespace GRINS
{

  class MeshManager
  {
    
  public:

    MeshManager();
    ~MeshManager();

    void read_input_options( const GetPot& input );

    // get/set pair for mesh object
    libMesh::Mesh* get_mesh(); 
    void set_mesh( libMesh::Mesh* mesh );

    void build_mesh();

  private:
    // all possible values for _mesh_option:
    //    "read_mesh_from_file", "create_1D_mesh", "create_2D_MESH", "create_3D_mesh", "mesh_already_loaded"
    std::string _mesh_option;
    bool _print_mesh_info_flag;
    std::string _mesh_filename;

    libMesh::Mesh* _mesh;

    // all possible values for _domain_type:
    //    "line", "rectangle", "cylinder", "box", sphere"
    std::string _domain_type;

    double _domain_x1_min, _domain_x2_min, _domain_x3_min;
    double _domain_x1_max, _domain_x2_max, _domain_x3_max;

    int _mesh_nx1, _mesh_nx2, _mesh_nx3;
    libMeshEnums::ElemType _element_type;

    bool _mesh_created_locally;
  };

} //End namespace block

#endif //GRINS_MESH_MANAGER_H
