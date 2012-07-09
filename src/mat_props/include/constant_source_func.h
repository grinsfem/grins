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

#ifndef CONSTANT_SOURCE_FUNC_H
#define CONSTANT_SOURCE_FUNC_H

// libMesh
#include "getpot.h"
#include "point.h"
#include "vector_value.h"

namespace GRINS
{
  class ConstantSourceFunction
  {
  public:

    ConstantSourceFunction( const GetPot& input );
    ~ConstantSourceFunction();

    inline
    Real operator()( libMesh::Point& ) const
    {
      return _value;
    }

    libMesh::Gradient grad( libMesh::Point& ) const
    {
      return libMesh::Gradient();
    }

  protected:

    Real _value;

  };
}
#endif //CONSTANT_SOURCE_FUNC_H
