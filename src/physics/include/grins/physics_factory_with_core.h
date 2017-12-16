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

#ifndef GRINS_PHYSICS_FACTORY_WITH_CORE_H
#define GRINS_PHYSICS_FACTORY_WITH_CORE_H

// GRINS
#include "physics_factory_base.h"

// C++
#include <map>

namespace GRINS
{
  //! PhysicsFactory base class for Physics that may have a related "core" Physics
  /*! There are some Physics subclasses that are inherently related to some "core"
    Physics. In such cases, some input options for the derived Physics class are
    slave to the "core" Physics. Thus, this provides a mechanism of naming the
    corresponding "core" Physics associated with the primary Physics. */
  class PhysicsFactoryWithCore : public PhysicsFactoryBase
  {
  public:

    PhysicsFactoryWithCore( const std::string& physics_name,
                            const std::string& core_physics_name );

    ~PhysicsFactoryWithCore(){};

  protected:

    std::string find_core_physics_name( const std::string& physics_name );

    //! Cache for "core" physics names
    /*! At parsing time, some Physics are slave to some input options from a core Physics.
      In such cases, we cache the core physics name. Note we use this function to
      avoid the static initialization fiasco. */
    static std::map<std::string,std::string>& core_physics_names();

  };

} // end namespace GRINS

#endif // GRINS_PHYSICS_FACTORY_WITH_CORE_H
