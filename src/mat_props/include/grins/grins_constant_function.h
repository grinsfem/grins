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

#ifndef GRINS_CONSTANT_FUNCTION_H
#define GRINS_CONSTANT_FUNCTION_H

// libMesh
#include "libmesh/libmesh.h"

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  class ConstantFunction
  {
  public:

    ConstantFunction();
    ~ConstantFunction();

    virtual void read_input_options( const GetPot& input ) =0;

    inline
    libMesh::Real operator()( libMesh::Real ) const
    { return _value; }

    inline
    libMesh::Real deriv( libMesh::Real ) const
    { return 0.0; }

  protected:

    libMesh::Real _value;
  };
}
#endif //GRINS_CONSTANT_FUNCTION_H
