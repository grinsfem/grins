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

#ifndef GRINS_PARSED_FUNCTION_DIRICHLET_OLD_STYLE_BC_FACTORY_H
#define GRINS_PARSED_FUNCTION_DIRICHLET_OLD_STYLE_BC_FACTORY_H

// GRINS
#include "grins/dirichlet_bc_factory_function_old_style_base.h"
#include "grins/parsed_function_traits.h"
#include "grins/parsed_function_factory_helper.h"

namespace GRINS
{
  template<typename FunctionType>
  class ParsedFunctionDirichletOldStyleBCFactory : public DirichletBCFactoryFunctionOldStyleBase<FunctionType>,
                                                   public ParsedFunctionFactoryHelper<FunctionType>
  {
  public:

    ParsedFunctionDirichletOldStyleBCFactory( const std::string& bc_type_name )
      : DirichletBCFactoryFunctionOldStyleBase<FunctionType>(bc_type_name),
      ParsedFunctionFactoryHelper<FunctionType>()
    {}

    ~ParsedFunctionDirichletOldStyleBCFactory(){};

  protected:

    virtual std::unique_ptr<FunctionType>
    build_func( const GetPot& input,
                MultiphysicsSystem& system,
                std::vector<std::string>& var_names,
                const std::string& section );

  };

  //! For notational convenience
  class ParsedDirichletOldStyleBCFactory : public ParsedFunctionDirichletOldStyleBCFactory<libMesh::FunctionBase<libMesh::Number> >
  {
  public:
    ParsedDirichletOldStyleBCFactory( const std::string& bc_type_name )
      : ParsedFunctionDirichletOldStyleBCFactory<libMesh::FunctionBase<libMesh::Number> >(bc_type_name)
    {}
  };

  //! For notational convenience
  class ParsedFEMDirichletOldStyleBCFactory : public ParsedFunctionDirichletOldStyleBCFactory<libMesh::FEMFunctionBase<libMesh::Number> >
  {
  public:
    ParsedFEMDirichletOldStyleBCFactory( const std::string& bc_type_name )
      : ParsedFunctionDirichletOldStyleBCFactory<libMesh::FEMFunctionBase<libMesh::Number> >(bc_type_name)
    {}
  };
} // end namespace GRINS

#endif // GRINS_PARSED_FUNCTION_DIRICHLET_OLD_STYLE_BC_FACTORY_H
