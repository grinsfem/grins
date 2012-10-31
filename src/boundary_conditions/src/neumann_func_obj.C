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

#include "neumann_func_obj.h"

GRINS::NeumannFuncObj::NeumannFuncObj( )
{
  return;
}

GRINS::NeumannFuncObj::~NeumannFuncObj( )
{
  return;
}

libMesh::Point GRINS::NeumannFuncObj::derivative( const libMesh::FEMContext& context, 
						  const unsigned int qp,
						  const GRINS::VariableIndex jac_var )
{
  // By default, does nothing.
  return libMesh::Point(0.0,0.0,0.0);
}

std::vector<GRINS::VariableIndex> GRINS::NeumannFuncObj:: get_other_jac_vars()
{
  return _jac_vars;
}
