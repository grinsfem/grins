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

#ifndef GRINS_PARSED_FUNCTION_DIRICHLET_BC_FACTORY_H
#define GRINS_PARSED_FUNCTION_DIRICHLET_BC_FACTORY_H

// GRINS
#include "grins/dirichlet_bc_factory_function_base.h"
#include "grins/parsed_function_traits.h"
#include "grins/parsed_function_factory_helper.h"

namespace GRINS
{
  template<typename FunctionType>
  class ParsedFunctionDirichletBCFactory : public DirichletBCFactoryFunctionBase<FunctionType>,
                                           public ParsedFunctionFactoryHelper<FunctionType>
  {
  public:

    ParsedFunctionDirichletBCFactory( const std::string& bc_type_name )
      : DirichletBCFactoryFunctionBase<FunctionType>(bc_type_name),
      ParsedFunctionFactoryHelper<FunctionType>()
    {}

    ~ParsedFunctionDirichletBCFactory(){};

  protected:
    //! Builds the Parsed(FEM)Function objects for boundary conditions
    /*! The variable names passed in will correspond to only a single
      VariableBase object, e.g. Velocity. The expected behavior is
      that if the user didn't specify a value for all the variables,
      then the unspecified variables will be set to zero. However,
      the user must've set at least one. The section arguments
      corresponds to the section to parse for the variables in the
      input file, e.g. input(section+"/"+var_names[0]). */
    virtual std::unique_ptr<FunctionType>
    build_func( const GetPot& input,
                MultiphysicsSystem& system,
                std::vector<std::string>& var_names,
                const std::string& section );

  };

  //! For notational convenience
  class ParsedDirichletBCFactory : public ParsedFunctionDirichletBCFactory<libMesh::FunctionBase<libMesh::Number> >
  {
  public:
    ParsedDirichletBCFactory( const std::string& bc_type_name )
      : ParsedFunctionDirichletBCFactory<libMesh::FunctionBase<libMesh::Number> >(bc_type_name)
    {}
  };

  //! For notational convenience
  class ParsedFEMDirichletBCFactory : public ParsedFunctionDirichletBCFactory<libMesh::FEMFunctionBase<libMesh::Number> >
  {
  public:
    ParsedFEMDirichletBCFactory( const std::string& bc_type_name )
      : ParsedFunctionDirichletBCFactory<libMesh::FEMFunctionBase<libMesh::Number> >(bc_type_name)
    {}
  };

} // end namespace GRINS

#endif // GRINS_PARSED_FUNCTION_DIRICHLET_BC_FACTORY_H
