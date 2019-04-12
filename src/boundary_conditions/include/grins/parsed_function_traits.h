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

#ifndef GRINS_PARSED_FUNCTION_TRAITS_H
#define GRINS_PARSED_FUNCTION_TRAITS_H

// libMesh
#include "libmesh/composite_fem_function.h"
#include "libmesh/composite_function.h"
#include "libmesh/const_fem_function.h"
#include "libmesh/fem_function_base.h"
#include "libmesh/function_base.h"
#include "libmesh/parsed_fem_function.h"
#include "libmesh/parsed_function.h"
#include "libmesh/zero_function.h"

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


  // Helper metafunctions
  template <typename FunctionType,
            bool is_fem_function =
            ParsedFunctionTraits<FunctionType>::is_fem_function>
  struct TypeFrom {
    typedef libMesh::CompositeFunction<libMesh::Number> to_composite;

    static libMesh::ParsedFunction<libMesh::Number>
    to_parsed(const libMesh::System & /* system */,
              const std::string & expression) {
      return libMesh::ParsedFunction<libMesh::Number>(expression);
    }

    static libMesh::ZeroFunction<libMesh::Number>
    to_zero() {
      return libMesh::ZeroFunction<libMesh::Number>();
    }
  };

  template <typename FunctionType>
  struct TypeFrom<FunctionType, true> {
    typedef libMesh::CompositeFEMFunction<libMesh::Number> to_composite;

    static libMesh::ParsedFEMFunction<libMesh::Number>
    to_parsed(const libMesh::System & system,
              const std::string & expression) {
      return libMesh::ParsedFEMFunction<libMesh::Number>(system, expression);
    }

    static libMesh::ConstFEMFunction<libMesh::Number>
    to_zero() {
      return libMesh::ConstFEMFunction<libMesh::Number>(0);
    }
  };

} // end namespace GRINS

#endif // GRINS_PARSED_FUNCTION_TRAITS_H
