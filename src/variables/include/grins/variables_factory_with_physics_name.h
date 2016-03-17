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

#ifndef GRINS_VARIABLES_FACTORY_WITH_PHYSICS_NAME_H
#define GRINS_VARIABLES_FACTORY_WITH_PHYSICS_NAME_H

// GRINS
#include "grins/variables_factory_base.h"

namespace GRINS
{
  //! Builds VariablesBase objects that require physics_name in the constructor
  /*! This factory will build VariablesBase subclasses that
      require a GetPot object to their constructor. Examples include
      the GenericVariable. */
  template<typename DerivedVariables>
  class VariablesFactoryWithPhysicsName : public VariablesFactoryBase
  {
  public:
    VariablesFactoryWithPhysicsName( const std::string& variable_name )
      : VariablesFactoryBase(variable_name)
    {}

    ~VariablesFactoryWithPhysicsName(){};

  protected:

    virtual libMesh::UniquePtr<VariablesBase> build_vars( const GetPot& input );

  };

} // end namespace GRINS

#endif // GRINS_VARIABLES_FACTORY_WITH_PHYSICS_NAME_H
