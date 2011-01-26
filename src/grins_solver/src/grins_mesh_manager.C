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

GRINS::MeshManager::MeshManager( const std::string mesh_options )
{
  _mesh_options = mesh_options;
  return;
}

GRINS::MeshManager::~MeshManager()
{
  return;
}

void GRINS::MeshManager::read_input_options( const GetPot& input )
{
  //TODO: get options
  //this->_mesh_input_option = input("mesh-options/mesh_input_option", 1 );

  return;
}

libMesh::Mesh* GRINS::MeshManager::get_mesh()
{
  return this->_mesh;
}

//TODO: discuss if this is needed
void GRINS::MeshManager::set_mesh( libMesh::Mesh* mesh )
{
  this->_mesh = mesh;
  return;
}
