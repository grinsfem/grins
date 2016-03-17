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

#ifndef GRINS_PHYSICS_FACTORY_INCOMPRESSIBLE_FLOW_H
#define GRINS_PHYSICS_FACTORY_INCOMPRESSIBLE_FLOW_H

// GRINS
#include "grins/physics_factory_with_core.h"

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

    virtual libMesh::UniquePtr<Physics> build_physics( const GetPot& input,
                                                       const std::string& physics_name );

    void visc_error_msg( const std::string& physics, const std::string& viscosity ) const;

  };

} // end namespace GRINS

#endif // GRINS_PHYSICS_FACTORY_INCOMPRESSIBLE_FLOW_H
