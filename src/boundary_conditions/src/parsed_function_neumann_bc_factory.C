//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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
#include "grins/parsed_function_neumann_bc_factory.h"

namespace GRINS
{
  ParsedFunctionNeumannBCFactory<libMesh::FunctionBase<libMesh::Number> >
  grins_factory_parsed_neumann("parsed_neumann");

  ParsedFunctionNeumannBCFactory<libMesh::FEMFunctionBase<libMesh::Number> >
  grins_factory_parsed_fem_neumann("parsed_fem_neumann");

  ParsedTractionBCFactory<libMesh::FunctionBase<libMesh::Number> >
  grins_factory_parsed_traction("parsed_traction");

  ParsedTractionBCFactory<libMesh::FEMFunctionBase<libMesh::Number> >
  grins_factory_parsed_fem_traction("parsed_fem_traction");

  // backward compatibility
  ParsedTractionBCFactory<libMesh::FunctionBase<libMesh::Number> >
  grins_factory_constant_traction("constant_traction");

} // end namespace GRINS
