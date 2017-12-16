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

#ifndef GRINS_PHYSICS_BUILDER_H
#define GRINS_PHYSICS_BUILDER_H

// C++
#include <string>

// GRINS
#include "grins/var_typedefs.h"

// C++
#include <set>

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  //! Manages runtime construction of Physics objects
  /*! This will parse the input file for the requested Physics
    and manage their construction. Actual construction of the
    Physics object is delegated to PhysicsFactoryBase subclasses.
    The PhysicsBuilder merely manages tasks around them as
    needed. To add a new Physics, the user should instantiate
    an appropriate PhysicsFactoryBase subclass. */
  class PhysicsBuilder
  {
  public:

    PhysicsBuilder(){};

    ~PhysicsBuilder(){};

    //! Returns container of constructed Physics objects
    static PhysicsList build_physics_map( const GetPot& input );

  private:

    //! Parses the requested Physics from input and populates the set passed to this function.
    static void parse_requested_physics( const GetPot& input, std::set<std::string>& requested_physics );

  };

} // end namespace GRINS

#endif // GRINS_PHYSICS_BUILDER_H
