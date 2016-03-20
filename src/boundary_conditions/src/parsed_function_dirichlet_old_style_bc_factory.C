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
  libMesh::UniquePtr<FunctionType>
  ParsedFunctionDirichletOldStyleBCFactory<FunctionType>::build_func( const GetPot& input,
                                                                      MultiphysicsSystem& system,
                                                                      std::vector<std::string>& var_names,
                                                                      const std::string& section )
  {
    libmesh_assert_equal_to( var_names.size(), 1 );

    std::vector<VariableIndex> dbc_vars(1,system.variable_number(var_names[0]));

    std::string section_str = section+"/"+DirichletBCFactoryFunctionOldStyleBase<FunctionType>::_value_var_old_style;
    std::string expression = input(section_str,"DIE!",DirichletBCFactoryFunctionOldStyleBase<FunctionType>::_value_idx_old_style);

    libMesh::UniquePtr<FunctionType> composite_func = this->build_composite_func();

    if( ParsedFunctionTraits<FunctionType>::is_fem_function )
      {
        libMesh::CompositeFEMFunction<libMesh::Number>* remapped_func =
          libMesh::libmesh_cast_ptr<libMesh::CompositeFEMFunction<libMesh::Number>*>(composite_func.get());

        libMesh::ParsedFEMFunction<libMesh::Number> parsed_func(system,expression);
        remapped_func->attach_subfunction(parsed_func, dbc_vars);
      }
    else
      {
        libMesh::CompositeFunction<libMesh::Number>* remapped_func =
          libMesh::libmesh_cast_ptr<libMesh::CompositeFunction<libMesh::Number>*>(composite_func.get());

        libMesh::ParsedFunction<libMesh::Number> parsed_func(expression);
        remapped_func->attach_subfunction(parsed_func, dbc_vars);
      }

    return composite_func;
  }

  // Instantiate all the ParsedDirichletOldStyle factories.
  ParsedDirichletOldStyleBCFactory grins_factory_parsed_dirichlet_old_style("parsed_dirichlet_old_style");
  ParsedDirichletOldStyleBCFactory grins_factory_constant_dirichlet_old_style("constant_dirichlet_old_style");
  ParsedFEMDirichletOldStyleBCFactory grins_factory_parsed_fem_dirichlet_old_style("parsed_fem_dirichlet_old_style");

} // end namespace GRINS
