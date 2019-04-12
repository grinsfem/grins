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
#include "grins/dirichlet_bc_factory_function_old_style_base.h"

namespace GRINS
{
  // Instantiate static members
  template<typename FunctionType>
  std::string DirichletBCFactoryFunctionOldStyleBase<FunctionType>::_value_var_old_style = std::string("DIE!");

  template<typename FunctionType>
  unsigned int DirichletBCFactoryFunctionOldStyleBase<FunctionType>::_value_idx_old_style = libMesh::invalid_uint;

  template<typename FunctionType>
  const std::vector<std::string>* DirichletBCFactoryFunctionOldStyleBase<FunctionType>::_var_names_old_style = NULL;

  template class GRINS::DirichletBCFactoryFunctionOldStyleBase<libMesh::FunctionBase<libMesh::Number> >;
  template class GRINS::DirichletBCFactoryFunctionOldStyleBase<libMesh::FEMFunctionBase<libMesh::Number> >;
} // end namespace GRINS
