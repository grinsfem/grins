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
// GRINS_Solver.h: Declarations for the GRINS_Solver class.
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_SOLVER_H
#define GRINS_SOLVER_H

#include "libmesh.h"

#include "GRINS_EquationSystems.h"

class GRINS_Solver
{

public:
  GRINS_Solver( const std::string application_options );
  ~GRINS_Solver();

private:
  std::string _application_options;
};

#endif
