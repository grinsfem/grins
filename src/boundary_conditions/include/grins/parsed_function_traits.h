//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_PARSED_FUNCTION_TRAITS_H
#define GRINS_PARSED_FUNCTION_TRAITS_H

// libMesh
#include "libmesh/function_base.h"
#include "libmesh/fem_function_base.h"

namespace GRINS
{
  template<typename FunctionType>
  struct ParsedFunctionTraits;

  template<typename FEShape>
  struct ParsedFunctionTraits<libMesh::FunctionBase<FEShape> >
  {
    static bool const is_fem_function = false;
  };

  template<typename FEShape>
  struct ParsedFunctionTraits<libMesh::FEMFunctionBase<FEShape> >
  {
    static bool const is_fem_function = true;
  };

} // end namespace GRINS

#endif // GRINS_PARSED_FUNCTION_TRAITS_H
