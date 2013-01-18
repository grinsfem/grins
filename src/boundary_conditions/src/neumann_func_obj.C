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

#include "grins/neumann_func_obj.h"

namespace GRINS
{

  NeumannFuncObj::NeumannFuncObj( )
  {
    return;
  }

  NeumannFuncObj::~NeumannFuncObj( )
  {
    return;
  }

  libMesh::Point NeumannFuncObj::derivative( const libMesh::FEMContext& context, 
					     const unsigned int qp,
					     const VariableIndex jac_var )
  {
    // By default, does nothing.
    return libMesh::Point(0.0,0.0,0.0);
  }

  std::vector<VariableIndex> NeumannFuncObj:: get_other_jac_vars()
  {
    return _jac_vars;
  }

} // namespace GRINS
