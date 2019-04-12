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

#ifndef GRINS_FACTORY_WITH_GETPOT_PHYSICS_NAME_H
#define GRINS_FACTORY_WITH_GETPOT_PHYSICS_NAME_H

// GRINS
#include "grins/factory_with_getpot.h"

namespace GRINS
{
  //! Abstract factory that provides availability of GetPot and a physics_name
  template<typename Base>
  class FactoryWithGetPotPhysicsName : public FactoryWithGetPot<Base>
  {
  public:
    FactoryWithGetPotPhysicsName( const std::string& name )
      : FactoryWithGetPot<Base>(name)
    {}

    ~FactoryWithGetPotPhysicsName(){};

    //! Setter for physics name
    /*! We need the physics_name to pass to the constructor, so we need
      to provide a hook to get it. Note that this should be the "full"
      physics name, including suffixes, etc. Subclasses dictate final
      behavior, but generally, this MUST be called each time build()
      is called as the expected behavior is for the physics_name to
      be reset after the build() call. */
    static void set_physics_name( const std::string& physics_name )
    { _physics_name = physics_name; }

  protected:

    static std::string _physics_name;

  };

} // end namespace GRINS

#endif // GRINS_FACTORY_WITH_GETPOT_PHYSICS_NAME_H
