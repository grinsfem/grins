//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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


#ifndef CONSTANT_SOURCE_FUNC_H
#define CONSTANT_SOURCE_FUNC_H

// GRINS
#include "grins/parameter_user.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/point.h"
#include "libmesh/vector_value.h"

namespace GRINS
{
  class ConstantSourceFunction : public ParameterUser
  {
  public:

    ConstantSourceFunction( const GetPot& input );

    ~ConstantSourceFunction() = default;

    inline
    libMesh::Real operator()( const libMesh::Point& ) const
    {
      return _value;
    }

    libMesh::Gradient grad( const libMesh::Point& ) const
    {
      return libMesh::Gradient();
    }

  protected:

    libMesh::Real _value;

  };
}
#endif //CONSTANT_SOURCE_FUNC_H
