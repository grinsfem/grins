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

#ifndef GRINS_NEUMANN_BC_PARSED_H
#define GRINS_NEUMANN_BC_PARSED_H

// GRINS
#include "grins/neumann_bc_function_base.h"
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/parsed_function.h"
#include "libmesh/composite_function.h"
#include "libmesh/parsed_fem_function.h"
#include "libmesh/composite_fem_function.h"

namespace GRINS
{
  template<typename FEShape>
  class ParsedNeumannBC : public NeumannBCFunctionBase<libMesh::FunctionBase<FEShape>,FEShape>
  {
  public:
    ParsedNeumannBC( const std::string& expression,
                     VariableIndex var )
      : NeumannBCFunctionBase<libMesh::FunctionBase<FEShape>,FEShape>( var )
    {
      this->_func.reset( new libMesh::ParsedFunction<FEShape>(expression) );
    }

    virtual ~ParsedNeumannBC(){};
  };

  template<typename FEShape>
  class CompositeParsedNeumannBC : public NeumannBCFunctionBase<libMesh::FunctionBase<FEShape>,FEShape>
  {
  public:
    CompositeParsedNeumannBC( const std::vector<std::string>& expressions,
                              const std::vector<VariableIndex>& vars )
      : NeumannBCFunctionBase<libMesh::FunctionBase<FEShape>,FEShape>( vars )
    {
      libmesh_assert_equal_to( expressions.size(), vars.size() );

      std::unique_ptr<libMesh::CompositeFunction<FEShape> >
        composite_func( new libMesh::CompositeFunction<FEShape> );

      for( unsigned int i = 0; i < vars.size(); i++ )
        {
          libMesh::ParsedFunction<FEShape> parsed_func(expressions[i]);
          std::vector<unsigned int> index(1,vars[i]);
          composite_func->attach_subfunction(parsed_func, index);
        }

      this->_func.reset(composite_func.release());
    }

    virtual ~CompositeParsedNeumannBC(){};
  };

  template<typename FEShape>
  class ParsedFEMNeumannBC : public NeumannBCFunctionBase<libMesh::FEMFunctionBase<FEShape>,FEShape>
  {
  public:
    ParsedFEMNeumannBC( const std::string& expression,
                        const MultiphysicsSystem& system,
                        VariableIndex var )
      : NeumannBCFunctionBase<libMesh::FEMFunctionBase<FEShape>,FEShape>( var )
    {
      this->_func.reset( new libMesh::ParsedFEMFunction<FEShape>(system,expression) );
    }

    virtual ~ParsedFEMNeumannBC(){};
  };

  template<typename FEShape>
  class CompositeParsedFEMNeumannBC : public NeumannBCFunctionBase<libMesh::FEMFunctionBase<FEShape>,FEShape>
  {
  public:
    CompositeParsedFEMNeumannBC( const std::vector<std::string>& expressions,
                                 const std::vector<VariableIndex>& vars,
                                 const MultiphysicsSystem& system )
      : NeumannBCFunctionBase<libMesh::FEMFunctionBase<FEShape>,FEShape>( vars )
    {
      libmesh_assert_equal_to( expressions.size(), vars.size() );

      std::unique_ptr<libMesh::CompositeFEMFunction<FEShape> >
        composite_func( new libMesh::CompositeFEMFunction<FEShape> );

      for( unsigned int i = 0; i < vars.size(); i++ )
        {
          libMesh::ParsedFEMFunction<FEShape> parsed_func(system,expressions[i]);
          std::vector<unsigned int> index(1,vars[i]);
          composite_func->attach_subfunction(parsed_func, index);
        }

      this->_func.reset(composite_func.release());
    }

    virtual ~CompositeParsedFEMNeumannBC(){};
  };
} // end namespace GRINS

#endif // GRINS_NEUMANN_BC_PARSED_H
