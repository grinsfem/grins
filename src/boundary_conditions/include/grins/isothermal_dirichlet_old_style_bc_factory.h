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

#ifndef GRINS_ISOTHERMAL_DIRICHLET_OLD_STYLE_BC_FACTORY_H
#define GRINS_ISOTHERMAL_DIRICHLET_OLD_STYLE_BC_FACTORY_H

// GRINS
#include "grins/dirichlet_bc_factory_function_old_style_base.h"

namespace GRINS
{
  class IsothermalDirichletOldStyleBCFactory : public DirichletBCFactoryFunctionOldStyleBase<libMesh::FunctionBase<libMesh::Number> >
  {
  public:

    IsothermalDirichletOldStyleBCFactory( const std::string& bc_type_name )
      : DirichletBCFactoryFunctionOldStyleBase<libMesh::FunctionBase<libMesh::Number> >(bc_type_name)
    {}

    ~IsothermalDirichletOldStyleBCFactory(){};

  protected:

    virtual libMesh::UniquePtr<libMesh::FunctionBase<libMesh::Number> >
    build_func( const GetPot& input,
                MultiphysicsSystem& system,
                std::vector<std::string>& var_names,
                const std::string& section );

  };

} // end namespace GRINS

#endif // GRINS_ISOTHERMAL_DIRICHLET_OLD_STYLE_BC_FACTORY_H
