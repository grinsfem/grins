//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
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
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef BC_TYPES_H
#define BC_TYPES_H

#include "dirichlet_func_obj.h"
#include "neumann_func_obj.h"

// Forward declare BC function objects.
//class GRINS::DirichletFuncObj;
//class GRINS::NeumannFuncObj;

namespace GRINS
{
  enum BC_TYPES{ NO_SLIP = 0,
		 PRESCRIBED_VELOCITY,
		 INFLOW,
		 AXISYMMETRIC,
		 ISOTHERMAL_WALL,
		 ADIABATIC_WALL,
		 PRESCRIBED_HEAT_FLUX,
		 GENERAL_HEAT_FLUX
  };

  typedef std::pair< GRINS::VariableIndex, std::tr1::shared_ptr<GRINS::DirichletFuncObj> > DBCMapPair;
  typedef std::map< GRINS::VariableIndex, std::tr1::shared_ptr<GRINS::DirichletFuncObj> > DirichletBCsMap;

  typedef std::pair< GRINS::VariableIndex, std::tr1::shared_ptr<GRINS::NeumannFuncObj> > NBCMapPair;
  typedef std::map< GRINS::VariableIndex, std::tr1::shared_ptr<GRINS::NeumannFuncObj> > NeumannBCsMap;

  typedef std::pair< GRINS::BoundaryID, GRINS::DirichletBCsMap > DBCContainerPair;
  typedef std::map< GRINS::BoundaryID, GRINS::DirichletBCsMap > DBCContainer;

  typedef std::pair< GRINS::BoundaryID, GRINS::NeumannBCsMap > NBCContainerPair;
  typedef std::map< GRINS::BoundaryID, GRINS::NeumannBCsMap > NBCContainer;
}
#endif //BC_TYPES_H
