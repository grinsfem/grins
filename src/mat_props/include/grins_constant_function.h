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

#ifndef GRINS_CONSTANT_FUNCTION_H
#define GRINS_CONSTANT_FUNCTION_H

// libMesh
#include "libmesh.h"
#include "getpot.h"

namespace GRINS
{
  class ConstantFunction
  {
  public:

    ConstantFunction();
    ~ConstantFunction();

    virtual void read_input_options( const GetPot& input ) =0;

    inline
    Real operator()( Real ) const
    { return _value; }

    inline
    Real deriv( Real ) const
    { return 0.0; }

  protected:

    Real _value;
  };
}
#endif //GRINS_CONSTANT_FUNCTION_H
