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

// This class
#include "grins/parsed_function_factory_helper.h"

// libMesh
#include "libmesh/composite_function.h"
#include "libmesh/composite_fem_function.h"
#include "libmesh/parsed_function.h"
#include "libmesh/parsed_fem_function.h"

namespace GRINS
{
  template<>
  std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >
  ParsedFunctionFactoryHelper<libMesh::FunctionBase<libMesh::Number> >::build_parsed_func
  ( const MultiphysicsSystem& /*system*/, const std::string& expression )
  {
    return std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >( new libMesh::ParsedFunction<libMesh::Number>(expression) );
  }

  template<>
  std::unique_ptr<libMesh::FEMFunctionBase<libMesh::Number> >
  ParsedFunctionFactoryHelper<libMesh::FEMFunctionBase<libMesh::Number> >::build_parsed_func
  ( const MultiphysicsSystem& system, const std::string& expression )
  {
    return std::unique_ptr<libMesh::FEMFunctionBase<libMesh::Number> >( new libMesh::ParsedFEMFunction<libMesh::Number>(system,expression) );
  }

  template<>
  std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >
  ParsedFunctionFactoryHelper<libMesh::FunctionBase<libMesh::Number> >::build_composite_func()
  {
    return std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >( new libMesh::CompositeFunction<libMesh::Number> );
  }

  template<>
  std::unique_ptr<libMesh::FEMFunctionBase<libMesh::Number> >
  ParsedFunctionFactoryHelper<libMesh::FEMFunctionBase<libMesh::Number> >::build_composite_func()
  {
    return std::unique_ptr<libMesh::FEMFunctionBase<libMesh::Number> >( new libMesh::CompositeFEMFunction<libMesh::Number> );
  }

} // end namespace GRINS
