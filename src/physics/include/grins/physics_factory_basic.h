//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_PHYSICS_FACTORY_BASIC_H
#define GRINS_PHYSICS_FACTORY_BASIC_H

#include "grins/physics_factory_base.h"

namespace GRINS
{
  template<typename DerivedPhysics>
  class PhysicsFactoryBasic : public PhysicsFactoryBase
  {
  public:
    PhysicsFactoryBasic( const std::string& physics_name )
      : PhysicsFactoryBase(physics_name)
    {}

    ~PhysicsFactoryBasic(){};

  protected:

    virtual std::unique_ptr<Physics> build_physics( const GetPot& input,
                                                    const std::string& physics_name )
    {
      return std::unique_ptr<Physics>( new DerivedPhysics(physics_name,input) );
    }

  };

} // end namespace GRINS

#endif // GRINS_PHYSICS_FACTORY_BASIC_H
