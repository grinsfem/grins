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
#ifndef POINT_FUNC_BASE_H
#define POINT_FUNC_BASE_H

#include "libmesh.h"
#include "point.h"

namespace GRINS
{

  //! Abstract base class for functions that return a libMesh::Point
  /*! Interface for all GRINS functions that return a point in space.
      Current applications are for Dirichlet boundary condition 
      applications. */
  class BasePointFuncObj
  {
  public:
    
    BasePointFuncObj();
    
    virtual ~BasePointFuncObj();
    
    virtual libMesh::Point operator()( const libMesh::Point& point ) = 0;
    
  }; // class BasePointFuncObj
  
} // namespace GRINS

#endif // POINT_FUNC_BASE_H
