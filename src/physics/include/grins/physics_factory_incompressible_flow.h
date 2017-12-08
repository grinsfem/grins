//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_PHYSICS_FACTORY_INCOMPRESSIBLE_FLOW_H
#define GRINS_PHYSICS_FACTORY_INCOMPRESSIBLE_FLOW_H

// GRINS
#include "grins/physics_factory_with_core.h"
#include "grins/physics_factory_helper.h"
#include "grins/constant_viscosity.h"
#include "grins/parsed_viscosity.h"
#include "grins/spalart_allmaras_viscosity.h"

namespace GRINS
{
  template<template<typename> class DerivedPhysics>
  class PhysicsFactoryIncompressibleFlow : public PhysicsFactoryWithCore
  {
  public:
    PhysicsFactoryIncompressibleFlow( const std::string& physics_name,
                                      const std::string& core_physics_name )
      : PhysicsFactoryWithCore(physics_name,core_physics_name)
    {}

    ~PhysicsFactoryIncompressibleFlow(){};

  protected:

    virtual std::unique_ptr<Physics> build_physics( const GetPot& input,
                                                       const std::string& physics_name );

    void visc_error_msg( const std::string& physics, const std::string& viscosity ) const;

  };

  template<template<typename> class DerivedPhysics>
  inline
  std::unique_ptr<Physics>
  PhysicsFactoryIncompressibleFlow<DerivedPhysics>::build_physics
  ( const GetPot& input, const std::string& physics_name )
  {
    std::string core_physics = this->find_core_physics_name(physics_name);

    std::string viscosity;
    PhysicsFactoryHelper::parse_viscosity_model(input,core_physics,viscosity);

    std::unique_ptr<Physics> new_physics;

    if( viscosity == "constant" )
      new_physics.reset( new DerivedPhysics<ConstantViscosity>(physics_name,input) );

    else if( viscosity == "parsed" )
      new_physics.reset( new DerivedPhysics<ParsedViscosity>(physics_name,input) );

    // For SA viscosity model, we need to parse what the "sub" viscosity model is
    else if( viscosity == "spalartallmaras" )
      {
        std::string turb_viscosity;
        PhysicsFactoryHelper::parse_turb_viscosity_model(input,core_physics,turb_viscosity);
        if( turb_viscosity == "constant" )
          new_physics.reset(new DerivedPhysics<SpalartAllmarasViscosity<ConstantViscosity> >(physics_name,input) );
        else
          this->visc_error_msg(physics_name, turb_viscosity);
      }
    else
      this->visc_error_msg(physics_name, viscosity);

    libmesh_assert(new_physics);

    return new_physics;
  }

  template<template<typename> class DerivedPhysics>
  inline
  void PhysicsFactoryIncompressibleFlow<DerivedPhysics>::visc_error_msg( const std::string& physics,
                                                                         const std::string& viscosity ) const
  {
    std::string error = "================================================================\n";
    error += "Invalid viscosity model for "+physics+"\n";
    error += "Viscosity model     = "+viscosity+"\n";
    error += "================================================================\n";

    libmesh_error_msg(error);
  }

} // end namespace GRINS

#endif // GRINS_PHYSICS_FACTORY_INCOMPRESSIBLE_FLOW_H
