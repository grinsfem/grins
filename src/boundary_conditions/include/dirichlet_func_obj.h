//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
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
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef DIRICHLET_FUNC_OBJ_H 
#define DIRICHLET_FUNC_OBJ_H

#include "libmesh.h"
#include "fem_context.h"

namespace GRINS
{

  //! Abstract base class for general, non-constant Dirichlet boundary conditions
  class DirichletFuncObj
  {
  public:
    
    DirichletFuncObj();
    
    virtual ~DirichletFuncObj();
    
    //! Returns the value of the implemented Dirichlet boundary condition
    /*! This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp. */
    virtual libMesh::Number value( const libMesh::FEMContext& context, const unsigned int qp ) = 0;

  }; // class DirichletFuncObj
  
} // namespace GRINS

#endif // DIRICHLET_FUNC_OBJ_H
