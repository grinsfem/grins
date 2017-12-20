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
#include "grins/physics_builder.h"

// GRINS
#include <memory>
#include "grins/physics_naming.h"
#include "grins/physics_factory_base.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  PhysicsList PhysicsBuilder::build_physics_map( const GetPot& input )
  {
    // First parse what Physics the user requested
    std::set<std::string> requested_physics;
    PhysicsBuilder::parse_requested_physics(input,requested_physics);

    // Give the PhysicsFactories access to the GetPot object
    PhysicsFactoryBase::set_getpot(input);

    PhysicsList physics_map;

    // Now populate the PhysicsList with with of the corresponding Physics objects
    for( std::set<std::string>::const_iterator physics_name_it = requested_physics.begin();
         physics_name_it != requested_physics.end();
         physics_name_it++ )
      {
        // Convenience
        std::string physics_name = *physics_name_it;

        // Set the physics_name in the PhysicsFactory as it will be used internally
        // to feed to the Physics constructors
        PhysicsFactoryBase::set_physics_name(physics_name);

        // Set the suffix in PhysicsNaming so it propagates to the Physics classes
        std::string physics_suffix = PhysicsNaming::extract_suffix(physics_name);
        PhysicsNaming::set_suffix(physics_suffix);

        // Strip any suffix out when querying the PhysicsFactories
        std::string physics_name_prefix = PhysicsNaming::extract_physics(physics_name);

        // Now build the actual Physics object
        std::unique_ptr<Physics> physics_ptr = PhysicsFactoryBase::build(physics_name_prefix);

        // Clear out the suffix now that we're done
        PhysicsNaming::clear_suffix();

        // Now give it to the PhysicsList
        physics_map[physics_name] = std::shared_ptr<Physics>( physics_ptr.release() );
      }

    // Echo Physics to display, if requested
    if( input( "screen-options/echo_physics", true ) )
      {
        libMesh::out << "==========================================================" << std::endl
                     << "List of Enabled Physics:" << std::endl;

        for( PhysicsListIter it = physics_map.begin();
             it != physics_map.end();
             it++ )
          {
            libMesh::out << it->first << std::endl;
          }
        libMesh::out <<  "==========================================================" << std::endl;
      }

    return physics_map;
  }

  void PhysicsBuilder::parse_requested_physics( const GetPot& input, std::set<std::string>& requested_physics )
  {
    int num_physics =  input.vector_variable_size("Physics/enabled_physics");

    if( num_physics < 1 )
      libmesh_error_msg("Error: Must enable at least one physics model!");

    for( int i = 0; i < num_physics; i++ )
      requested_physics.insert( input("Physics/enabled_physics", "NULL", i ) );
  }

} // end namespace GRINS
