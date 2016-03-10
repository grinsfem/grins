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

#ifndef GRINS_PHYSICS_FACTORY_BASE_H
#define GRINS_PHYSICS_FACTORY_BASE_H

// GRINS
#include "grins/factory_with_getpot.h"
#include "grins/physics.h"

namespace GRINS
{
  //! Builds Physics objects, used by PhysicsBuilder
  /*! The user may subclass this class for more building more complex Physics objects.
      Because Physics objects require a GetPot input file object and the physics_name
      at construction time, both set_getpot() and  set_physics_name() MUST be called
      before build() function. Note that set_physics_name() MUST be called each time
      a new Physics is built.*/
  class PhysicsFactoryBase : public FactoryWithGetPot<Physics>
  {
  public:
    PhysicsFactoryBase( const std::string& physics_name )
      : FactoryWithGetPot<Physics>(physics_name)
    {}

    ~PhysicsFactoryBase(){};

    //! Setter for physics name
    /*! We need the physics_name to pass to the constructor, so we need
        to provide a hook to get it. Note that this should be the "full"
        physics name, including suffixes, etc.  MUST be called each time
        a new Physics is built as we reset the internal name each time. */
    static void set_physics_name( const std::string& physics_name )
    { _physics_name = physics_name; }

  protected:

    virtual libMesh::UniquePtr<Physics> build_physics( const GetPot& input,
                                                       const std::string& physics_name ) =0;

    static std::string _physics_name;

  private:

    virtual libMesh::UniquePtr<Physics> create();

  };

  inline
  libMesh::UniquePtr<Physics> PhysicsFactoryBase::create()
  {
    // Make sure user set the physics name
    if( _physics_name == std::string("DIE!") )
      libmesh_error_msg("ERROR: must call set_physics_name() before building Physics!");

    if( !_input )
      libmesh_error_msg("ERROR: must call set_getpot() before building Physics!");

    libMesh::UniquePtr<Physics> new_physics = this->build_physics( *_input, _physics_name );

    // Reset the _physics_name for error checking
    _physics_name = std::string("DIE!");

    libmesh_assert(new_physics);

    return new_physics;
  }

} // end namespace GRINS

#endif // GRINS_PHYSICS_FACTORY_BASE_H
