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

// This class
#include "grins/physics_factory_incompressible_turb_flow.h"

// GRINS
#include "grins/physics_factory_helper.h"
#include "grins/constant_viscosity.h"

// Physics whose factories we're instantiating
#include "grins/spalart_allmaras.h"
#include "grins/spalart_allmaras_spgsm_stab.h"

namespace GRINS
{
  template<template<typename> class DerivedPhysics>
  libMesh::UniquePtr<Physics> PhysicsFactoryIncompressibleTurbFlow<DerivedPhysics>::build_physics( const GetPot& input,
                                                                                                   const std::string& physics_name )
  {
    std::string core_physics = this->find_core_physics_name(physics_name);

    std::string viscosity;
    PhysicsFactoryHelper::parse_turb_viscosity_model(input,core_physics,viscosity);

    libMesh::UniquePtr<Physics> new_physics;

    if( viscosity == "constant" )
      new_physics.reset( new DerivedPhysics<ConstantViscosity>(physics_name,input) );

    else
      this->visc_error_msg(physics_name, viscosity);

    libmesh_assert(new_physics);

    return new_physics;
  }

  template<template<typename> class DerivedPhysics>
  void PhysicsFactoryIncompressibleTurbFlow<DerivedPhysics>::visc_error_msg( const std::string& physics,
                                                                             const std::string& viscosity ) const
  {
    std::string error = "================================================================\n";
    error += "Invalid turblence viscosity model for "+physics+"\n";
    error += "Viscosity model     = "+viscosity+"\n";
    error += "================================================================\n";

    libmesh_error_msg(error);
  }

  // Instantiate all the "turbulent incompressble flow" Physics factories.
  PhysicsFactoryIncompressibleTurbFlow<SpalartAllmaras> grins_factory_spalart_allmaras
  (PhysicsNaming::spalart_allmaras(),PhysicsNaming::spalart_allmaras());

  PhysicsFactoryIncompressibleTurbFlow<SpalartAllmarasSPGSMStabilization> grins_factory_spalart_allmaras_spgsm_stab
  (PhysicsNaming::spalart_allmaras_spgsm_stab(),PhysicsNaming::spalart_allmaras());

} // end namespace GRINS
