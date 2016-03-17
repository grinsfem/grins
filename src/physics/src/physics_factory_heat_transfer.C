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
#include "grins/physics_factory_heat_transfer.h"

// GRINS
#include "grins/physics_factory_helper.h"
#include "grins/constant_conductivity.h"
#include "grins/parsed_conductivity.h"

// Physics we're instantiating
#include "grins/heat_conduction.h"
#include "grins/heat_transfer.h"
#include "grins/heat_transfer_adjoint_stab.h"
#include "grins/heat_transfer_spgsm_stab.h"
#include "grins/axisym_heat_transfer.h"

namespace GRINS
{
  template<template<typename> class DerivedPhysics>
  libMesh::UniquePtr<Physics> PhysicsFactoryHeatTransfer<DerivedPhysics>::build_physics( const GetPot& input,
                                                                                         const std::string& physics_name )
  {
    std::string core_physics = this->find_core_physics_name(physics_name);

    std::string conductivity;
    PhysicsFactoryHelper::parse_conductivity_model(input,core_physics,conductivity);

    libMesh::UniquePtr<Physics> new_physics;

    if( conductivity == "constant" )
      new_physics.reset( new DerivedPhysics<ConstantConductivity>(physics_name,input) );

    else if( conductivity == "parsed" )
      new_physics.reset( new DerivedPhysics<ParsedConductivity>(physics_name,input) );

    else
      this->cond_error_msg(physics_name, conductivity);

    libmesh_assert(new_physics);

    return new_physics;
  }

  template<template<typename> class DerivedPhysics>
  void PhysicsFactoryHeatTransfer<DerivedPhysics>::cond_error_msg( const std::string& physics,
                                                                   const std::string& conductivity ) const
  {
    std::string error = "================================================================\n";
    error += "Invalid conductivity model for "+physics+"\n";
    error += "Conductivity model     = "+conductivity+"\n";
    error += "================================================================\n";

    libmesh_error_msg(error);
  }

  PhysicsFactoryHeatTransfer<HeatConduction> grins_factory_heat_conduction
  (PhysicsNaming::heat_conduction(),PhysicsNaming::heat_conduction());

  PhysicsFactoryHeatTransfer<HeatTransfer> grins_factory_heat_transfer
  (PhysicsNaming::heat_transfer(),PhysicsNaming::heat_transfer());

  PhysicsFactoryHeatTransfer<HeatTransferAdjointStabilization> grins_factory_heat_transfer_adjoint_stab
  (PhysicsNaming::heat_transfer_adjoint_stab(),PhysicsNaming::heat_transfer());

  PhysicsFactoryHeatTransfer<HeatTransferAdjointStabilization> grins_factory_heat_transfer_spgsm_stab
  (PhysicsNaming::heat_transfer_spgsm_stab(),PhysicsNaming::heat_transfer());

  // This needs to die. Axisymmetry should be handled within heat_transfer
  PhysicsFactoryHeatTransfer<AxisymmetricHeatTransfer> grins_factory_axi_heat_transfer
  (PhysicsNaming::axisymmetric_heat_transfer(),PhysicsNaming::heat_transfer());

} // end namespace GRINS
