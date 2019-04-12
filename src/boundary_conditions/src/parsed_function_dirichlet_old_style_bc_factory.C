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
#include "grins/parsed_function_dirichlet_old_style_bc_factory.h"

// libMesh
#include "libmesh/composite_function.h"
#include "libmesh/composite_fem_function.h"
#include "libmesh/parsed_function.h"
#include "libmesh/parsed_fem_function.h"

namespace GRINS
{
  template<typename FunctionType>
  std::unique_ptr<FunctionType>
  ParsedFunctionDirichletOldStyleBCFactory<FunctionType>::build_func( const GetPot& input,
                                                                      MultiphysicsSystem& system,
                                                                      std::vector<std::string>& var_names,
                                                                      const std::string& section )
  {
    libmesh_assert_equal_to( var_names.size(), 1 );

    std::vector<VariableIndex> dbc_vars(1,system.variable_number(var_names[0]));

    std::string section_str = section+"/"+DirichletBCFactoryFunctionOldStyleBase<FunctionType>::_value_var_old_style;
    std::string expression = input(section_str,"DIE!",DirichletBCFactoryFunctionOldStyleBase<FunctionType>::_value_idx_old_style);

    std::unique_ptr<FunctionType> all_funcs = this->build_composite_func();

    typedef typename TypeFrom<FunctionType>::to_composite composite_type;
    composite_type * composite_func =
      libMesh::cast_ptr<composite_type *>(all_funcs.get());

    composite_func->attach_subfunction
      (TypeFrom<FunctionType>::to_parsed(system, expression), dbc_vars);

    return all_funcs;
  }

  // Instantiate
  template class ParsedFunctionDirichletOldStyleBCFactory<libMesh::FunctionBase<libMesh::Number> >;
  template class ParsedFunctionDirichletOldStyleBCFactory<libMesh::FEMFunctionBase<libMesh::Number> >;

} // end namespace GRINS
