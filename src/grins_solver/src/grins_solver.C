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
// Definitions for the GRINS::Solver class.
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "grins_solver.h"

#include <iostream>

GRINS::Solver::Solver( const std::string application_options )
  : _output_vis_flag(false)
{
  std::cout << " GRINS::Solver constructor ..." << std::endl;
  _application_options = application_options;
  return;
}

GRINS::Solver::~Solver()
{
  std::cout << " GRINS::Solver  destructor ..." << std::endl;
  return;
}

void GRINS::Solver::read_input_options( const GetPot& input )
{
  this->_output_vis_flag = input("vis-options/output_vis_flag", false );

  //TODO: Currently there for quick and stupid test. Delete.
  std::cout << "_output_vis_flag value = " << _output_vis_flag << std::endl;

  return;
}

libMesh::Mesh* GRINS::Solver::get_mesh()
{
  return this->_mesh;
}

void GRINS::Solver::set_mesh( libMesh::Mesh *mesh )
{
  this->_mesh = mesh;
  return;
}
