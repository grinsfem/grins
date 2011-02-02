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
// Definitions for the GRINS::MeshManager class.
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "grins_mesh_manager.h"

#include <iostream>

GRINS::MeshManager::MeshManager()
  : _mesh(NULL)
{
  return;
}

GRINS::MeshManager::~MeshManager()
{
  return;
}

void GRINS::MeshManager::read_input_options( const GetPot& input )
{
  // TODO: should be fine to have (int)READ_MESH_FROM_FILE
  this->_mesh_option = input("mesh-options/mesh_option", (int)READ_MESH_FROM_FILE );
  this->_print_mesh_info_flag = input("mesh-options/print_mesh_info_flag", false );

  return;
}

libMesh::Mesh* GRINS::MeshManager::get_mesh()
{
  return this->_mesh;
}

void GRINS::MeshManager::set_mesh( libMesh::Mesh* mesh )
{
  this->_mesh = mesh;
  return;
}

void GRINS::MeshManager::build_mesh()
{

  if( this->_mesh_option==MESH_ALREADY_LOADED )
    {
      std::cerr << " GRINS::MeshManager::build_mesh() : mesh already loaded " << std::endl;
      exit(1); // TODO: something more sophisticated for parallel runs?
    }

  switch (this->_mesh_option)
  {
    case READ_MESH_FROM_FILE:
      {
        // TODO: fill
      }
      break;
    case CREATE_1D_MESH:
      {
        // TODO: fill
      }
      break;
    case CREATE_2D_MESH:
      {
        // TODO: fill
      }
      break;
    case CREATE_3D_MESH:
      {
        // TODO: fill
      }
      break;
    default:
      {
        // TODO: fill
      }
      break;
  }

  return;
}
