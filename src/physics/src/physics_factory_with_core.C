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

// This class
#include "grins/physics_factory_with_core.h"

namespace GRINS
{
  PhysicsFactoryWithCore::PhysicsFactoryWithCore( const std::string& physics_name,
                                                  const std::string& core_physics_name )
    : PhysicsFactoryBase(physics_name)
  {
    this->core_physics_names().insert( std::make_pair(physics_name, core_physics_name) );
  }

  std::string PhysicsFactoryWithCore::find_core_physics_name( const std::string& physics_name )
  {
    std::string physics_name_stripped = PhysicsNaming::extract_physics(physics_name);

    if( this->core_physics_names().find(physics_name_stripped) == this->core_physics_names().end() )
      libmesh_error_msg("ERROR: Could not find core_physics associated with "+physics_name_stripped+"!");

    std::string core_physics_no_suffix = this->core_physics_names().find(physics_name_stripped)->second;

    // The user will be expecting any suffix to be on the returned core physics name
    std::string core_physics = core_physics_no_suffix+PhysicsNaming::extract_suffix(physics_name);

    return core_physics;
  }

  std::map<std::string,std::string>& PhysicsFactoryWithCore::core_physics_names()
  {
    static std::map<std::string,std::string> _core_physics_names;
    return _core_physics_names;
  }

} // end namespace GRINS
