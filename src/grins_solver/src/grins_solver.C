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
// Definitions for the GRINS_Solver class.
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "grins_solver.h"

#include <iostream>

GRINS::GRINSSolver::GRINSSolver( const std::string application_options )
{
  std::cout << " GRINS_Solver constructor ..." << std::endl;
  _application_options = application_options;
  return;
}

GRINS::GRINSSolver::~GRINSSolver()
{
  std::cout << " GRINS_Solver  destructor ..." << std::endl;
  return;
}
