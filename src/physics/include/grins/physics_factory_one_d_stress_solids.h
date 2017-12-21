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

#ifndef GRINS_PHYSICS_FACTORY_ONE_D_STRESS_SOLIDS_H
#define GRINS_PHYSICS_FACTORY_ONE_D_STRESS_SOLIDS_H

// GRINS
#include "grins/materials_parsing.h"
#include "grins/physics_factory_with_core.h"
#include "grins/hookes_law_1d.h"

namespace GRINS
{
  template<template<typename> class DerivedPhysics>
  class PhysicsFactoryOneDStressSolids : public PhysicsFactoryWithCore
  {
  public:
    PhysicsFactoryOneDStressSolids( const std::string& physics_name,
                                    const std::string& core_physics_name )
      : PhysicsFactoryWithCore(physics_name,core_physics_name)
    {}

    ~PhysicsFactoryOneDStressSolids(){};

  protected:

    virtual std::unique_ptr<Physics> build_physics( const GetPot& input,
                                                    const std::string& physics_name );

  };

  template<template<typename> class DerivedPhysics>
  inline
  std::unique_ptr<Physics> PhysicsFactoryOneDStressSolids<DerivedPhysics>::build_physics( const GetPot& input,
                                                                                          const std::string& physics_name )
  {
    std::string core_physics = this->find_core_physics_name(physics_name);

    std::string model = "none";
    std::string strain_energy = "none";

    MaterialsParsing::stress_strain_model( input, core_physics, model, strain_energy );

    std::unique_ptr<Physics> new_physics;

    if( model == std::string("hookes_law") )
      {
        new_physics.reset( new DerivedPhysics<HookesLaw1D>
                           (physics_name,input, false /*is_compressible*/) );
      }
    else
      {
        std::string error = "Error: Invalid stress-strain model: "+model+"!\n";
        error += "       Valid values are: hookes_law\n";
        libmesh_error_msg(error);
      }

    libmesh_assert(new_physics);

    return new_physics;
  }

} // end namespace GRINS

#endif // GRINS_PHYSICS_FACTORY_ONE_D_STRESS_SOLIDS_H
