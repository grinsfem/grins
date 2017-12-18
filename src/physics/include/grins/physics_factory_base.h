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

#ifndef GRINS_PHYSICS_FACTORY_BASE_H
#define GRINS_PHYSICS_FACTORY_BASE_H

// GRINS
#include "grins/factory_with_getpot_physics_name.h"
#include "grins/physics.h"

namespace GRINS
{
  // According to the standard, we need a declaration of the
  // specialization which precedes any automatic instantiation.
  template<> std::string FactoryWithGetPotPhysicsName<Physics>::_physics_name;
  template<> const GetPot* FactoryWithGetPot<Physics>::_input;

  //! Builds Physics objects, used by PhysicsBuilder
  /*! The user may subclass this class for more building more complex Physics objects.
    Because Physics objects require a GetPot input file object and the physics_name
    at construction time, both set_getpot() and  set_physics_name() MUST be called
    before build() function. Note that set_physics_name() MUST be called each time
    a new Physics is built.*/
  class PhysicsFactoryBase : public FactoryWithGetPotPhysicsName<Physics>
  {
  public:
    PhysicsFactoryBase( const std::string& physics_name )
      : FactoryWithGetPotPhysicsName<Physics>(physics_name)
    {}

    ~PhysicsFactoryBase(){};

  protected:

    virtual std::unique_ptr<Physics> build_physics( const GetPot& input,
                                                    const std::string& physics_name ) =0;

  private:

    virtual std::unique_ptr<Physics> create();

  };

  inline
  std::unique_ptr<Physics> PhysicsFactoryBase::create()
  {
    // Make sure user set the physics name
    if( _physics_name == std::string("DIE!") )
      libmesh_error_msg("ERROR: must call set_physics_name() before building Physics!");

    if( !_input )
      libmesh_error_msg("ERROR: must call set_getpot() before building Physics!");

    std::unique_ptr<Physics> new_physics = this->build_physics( *_input, _physics_name );

    // Reset the _physics_name for error checking
    _physics_name = std::string("DIE!");

    libmesh_assert(new_physics);

    return new_physics;
  }

} // end namespace GRINS

#endif // GRINS_PHYSICS_FACTORY_BASE_H
