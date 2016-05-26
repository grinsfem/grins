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
#include "grins/parsed_function_dirichlet_bc_factory.h"

// libMesh
#include "libmesh/composite_function.h"
#include "libmesh/composite_fem_function.h"
#include "libmesh/parsed_function.h"
#include "libmesh/parsed_fem_function.h"
#include "libmesh/zero_function.h"
#include "libmesh/const_fem_function.h"

namespace GRINS
{
  template<typename FunctionType>
  libMesh::UniquePtr<FunctionType>
  ParsedFunctionDirichletBCFactory<FunctionType>::build_func( const GetPot& input,
                                                              MultiphysicsSystem& system,
                                                              std::vector<std::string>& var_names,
                                                              const std::string& section )
  {
    libmesh_assert( !var_names.empty() );

    //! This is really a "composite" function. We'll cast in the helper functions.
    libMesh::UniquePtr<FunctionType> all_funcs;

    // We're given the active variables in var_names. Let's first check
    // which ones the user actually set in the input file.
    // If there's only one variable in var_names, then check_for_vars will
    // error if the user didn't set it, so we don't need to consider that
    // case here.
    std::set<std::string> vars_found;
    this->check_for_vars(input,section,var_names,&vars_found);

    all_funcs.reset( this->build_composite_func().release() );

    for( std::vector<std::string>::const_iterator var = var_names.begin();
         var < var_names.end(); ++var )
      {
        std::vector<VariableIndex> var_idx(1,system.variable_number(*var));

        // If the user set this variable in input, parse the expression and
        // add it to the composition function
        if( vars_found.find(*var) != vars_found.end() )
          {
            std::string expression = input(section+"/"+(*var),"DIE!");

            if( ParsedFunctionTraits<FunctionType>::is_fem_function )
              {
                libMesh::CompositeFEMFunction<libMesh::Number>* composite_func =
                  libMesh::libmesh_cast_ptr<libMesh::CompositeFEMFunction<libMesh::Number>*>(all_funcs.get());

                libMesh::ParsedFEMFunction<libMesh::Number> parsed_func(system,expression);
                composite_func->attach_subfunction(parsed_func,var_idx);
              }
            else
              {
                libMesh::CompositeFunction<libMesh::Number>* composite_func =
                  libMesh::libmesh_cast_ptr<libMesh::CompositeFunction<libMesh::Number>*>(all_funcs.get());

                libMesh::ParsedFunction<libMesh::Number> parsed_func(expression);
                composite_func->attach_subfunction(parsed_func,var_idx);
              }
          }
        // Otherwise, we set this variable to be zero.
        else
          {
            if( ParsedFunctionTraits<FunctionType>::is_fem_function )
              {
                libMesh::CompositeFEMFunction<libMesh::Number>* composite_func =
                  libMesh::libmesh_cast_ptr<libMesh::CompositeFEMFunction<libMesh::Number>*>(all_funcs.get());

                libMesh::ConstFEMFunction<libMesh::Number> zero_func(0.0);
                composite_func->attach_subfunction(zero_func,var_idx);
              }
            else
              {
                libMesh::CompositeFunction<libMesh::Number>* composite_func =
                  libMesh::libmesh_cast_ptr<libMesh::CompositeFunction<libMesh::Number>*>(all_funcs.get());

                libMesh::ZeroFunction<libMesh::Number> zero_func;
                composite_func->attach_subfunction(zero_func,var_idx);
              }
          }
      }

    return all_funcs;
  }

  // Instantiate all the Parsed(FEM)Function factories.
  ParsedDirichletBCFactory grins_factory_parsed_dirichlet("parsed_dirichlet");
  ParsedFEMDirichletBCFactory grins_factory_parsed_fem_dirichlet("parsed_fem_dirichlet");
  ParsedDirichletBCFactory grins_factory_parsed_displacement("parsed_displacement");

} // end namespace GRINS
