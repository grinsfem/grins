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

#ifndef GRINS_HOMOGENEOUS_DIRICHLET_BC_FACTORY_H
#define GRINS_HOMOGENEOUS_DIRICHLET_BC_FACTORY_H

// GRINS
#include "grins/dirichlet_bc_factory_function_base.h"

// libMesh
#include "libmesh/zero_function.h"

namespace GRINS
{
  class HomogeneousDirichletBCFactory : public DirichletBCFactoryFunctionBase<libMesh::FunctionBase<libMesh::Number> >
  {
  public:

    HomogeneousDirichletBCFactory( const std::string& bc_type_name )
      : DirichletBCFactoryFunctionBase(bc_type_name)
    {}

    ~HomogeneousDirichletBCFactory(){};

  protected:

    //! All the variables are 0, so just return 0 function.
    virtual std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >
    build_func( const GetPot& /*input*/,
                MultiphysicsSystem& /*system*/,
                std::vector<std::string>& /*var_names*/,
                const std::string& /*section*/ )
    {
      return std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >( new libMesh::ZeroFunction<libMesh::Number> );
    }

  };

} // end namespace GRINS

#endif // GRINS_HOMOGENEOUS_DIRICHLET_BC_FACTORY_H
