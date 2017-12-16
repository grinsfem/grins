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
#include "grins/physics_naming.h"

namespace GRINS
{

  void PhysicsNaming::extract_physics_and_suffix( const std::string& full_name,
                                                  std::string& physics_name,
                                                  std::string& suffix )
  {
    physics_name = PhysicsNaming::extract_physics(full_name);
    suffix = PhysicsNaming::extract_suffix(full_name);
  }

  std::string PhysicsNaming::extract_physics( const std::string& full_name )
  {
    // We look for this delimiter to separate the physics name from
    // the user-supplied suffix
    std::string delimiter = PhysicsNaming::physics_name_delimiter();

    std::size_t idx = full_name.find_first_of( delimiter );

    std::string physics_name;

    // If this delimiter is not found, then the full_name is the physics_name
    if( idx == full_name.npos )
      physics_name = full_name;

    // If it is found, the first part is the physics_name
    // and the second part is the suffix
    else
      physics_name = full_name.substr(0,idx);

    return physics_name;
  }

  std::string PhysicsNaming::extract_suffix( const std::string& full_name )
  {
    // We look for this delimiter to separate the physics name from
    // the user-supplied suffix
    std::string delimiter = PhysicsNaming::physics_name_delimiter();

    std::size_t idx = full_name.find_first_of( delimiter );

    std::string suffix;

    // If this delimiter is not found, then there's no suffix
    // So no need to do anything in that case.
    // If it is found, the first part is the physics_name
    // and the second part is the suffix.
    // Note the suffix extraction *includes* the delimiter
    if( idx != full_name.npos )
      suffix = full_name.substr(idx,full_name.npos);

    return suffix;
  }

} // end namespace GRINS
