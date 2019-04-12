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
#include "grins/parsed_function_neumann_bc_factory_helper.h"

// GRINS
#include "grins/neumann_bc_parsed.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  template<typename FunctionType>
  std::shared_ptr<NeumannBCAbstract>
  ParsedFunctionNeumannBCFactoryHelper<FunctionType>::build_neumman_func_common( const GetPot& input,
                                                                                 MultiphysicsSystem& system,
                                                                                 const FEVariablesBase& fe_var,
                                                                                 const std::string& flux_input )
  {
    const std::vector<std::string>& var_names = fe_var.active_var_names();

    std::shared_ptr<NeumannBCAbstract> func;

    // Use "standard" parsed version if there's only one variable
    if( var_names.size() == 1 )
      {
        libmesh_assert_equal_to( fe_var.var_indices().size(), 1 );
        std::string expression = input(flux_input,std::string("DIE!"));
        func = this->build_parsed_neumann_func(system,expression,fe_var.var_indices()[0]);
      }
    // Otherwise, use the composite versions
    else
      {
        libmesh_assert_equal_to( fe_var.var_indices().size(), var_names.size() );

        // We already checked size consistency for flux input and var_names
        // so just use var_names for the size
        std::vector<std::string> expressions(var_names.size());

        for( unsigned int i = 0; i < var_names.size(); i++ )
          expressions[i] = input(flux_input,std::string("DIE!"),i);

        func = this->build_composite_parsed_neumann_func(system,expressions,fe_var.var_indices());
      }

    return func;
  }

  template<>
  std::shared_ptr<NeumannBCAbstract>
  ParsedFunctionNeumannBCFactoryHelper<libMesh::FunctionBase<libMesh::Number> >::build_parsed_neumann_func
  (MultiphysicsSystem& /*system*/, const std::string& expression, VariableIndex var_idx )
  {
    return std::shared_ptr<NeumannBCAbstract>( new ParsedNeumannBC<libMesh::Number>(expression,var_idx) );
  }

  template<>
  std::shared_ptr<NeumannBCAbstract>
  ParsedFunctionNeumannBCFactoryHelper<libMesh::FEMFunctionBase<libMesh::Number> >::build_parsed_neumann_func
  (MultiphysicsSystem& system, const std::string& expression, VariableIndex var_idx )
  {
    return std::shared_ptr<NeumannBCAbstract>( new ParsedFEMNeumannBC<libMesh::Number>(expression,system,var_idx) );
  }

  template<>
  std::shared_ptr<NeumannBCAbstract>
  ParsedFunctionNeumannBCFactoryHelper<libMesh::FunctionBase<libMesh::Number> >::build_composite_parsed_neumann_func
  (MultiphysicsSystem& /*system*/, const std::vector<std::string>& expressions,
   const std::vector<VariableIndex>& var_indices )
  {
    return std::shared_ptr<NeumannBCAbstract>( new CompositeParsedNeumannBC<libMesh::Number>(expressions,var_indices) );
  }

  template<>
  std::shared_ptr<NeumannBCAbstract>
  ParsedFunctionNeumannBCFactoryHelper<libMesh::FEMFunctionBase<libMesh::Number> >::build_composite_parsed_neumann_func
  (MultiphysicsSystem& system, const std::vector<std::string>& expressions,
   const std::vector<VariableIndex>& var_indices )
  {
    return std::shared_ptr<NeumannBCAbstract>( new CompositeParsedFEMNeumannBC<libMesh::Number>(expressions,var_indices,system) );
  }

  // Instantiate
  template class ParsedFunctionNeumannBCFactoryHelper<libMesh::FunctionBase<libMesh::Number> >;
  template class ParsedFunctionNeumannBCFactoryHelper<libMesh::FEMFunctionBase<libMesh::Number> >;

} // end namespace GRINS
