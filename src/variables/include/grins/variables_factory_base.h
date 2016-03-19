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

#ifndef GRINS_VARIABLES_FACTORY_BASE_H
#define GRINS_VARIABLES_FACTORY_BASE_H

// GRINS
#include "grins/factory_with_getpot_physics_name.h"
#include "grins/variables_base.h"

namespace GRINS
{
  //! Builds VariableBase objects
  /*! Most variable classes only require a GetPot object to the constructor, but
      others may require more information. Subclasses can dictate the necessary
      behavior. */
  class VariablesFactoryBase : public FactoryWithGetPotPhysicsName<VariablesBase>
  {
  public:
    VariablesFactoryBase( const std::string& variable_name )
      : FactoryWithGetPotPhysicsName<VariablesBase>(variable_name)
    {}

    ~VariablesFactoryBase(){};

  protected:

    //! Subclasses implement this method for building the VariablesBase object.
    virtual libMesh::UniquePtr<VariablesBase> build_vars( const GetPot& input ) =0;

  private:

    virtual libMesh::UniquePtr<VariablesBase> create();

  };

  inline
  libMesh::UniquePtr<VariablesBase> VariablesFactoryBase::create()
  {
    if( !_input )
      libmesh_error_msg("ERROR: must call set_getpot() before building VariablesBase!");

    libMesh::UniquePtr<VariablesBase> new_vars = this->build_vars( *_input );

    libmesh_assert(new_vars);

    return new_vars;
  }

} // end namespace GRINS

#endif // GRINS_VARIABLES_FACTORY_BASE_H
