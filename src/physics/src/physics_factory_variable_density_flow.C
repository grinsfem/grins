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

// This class
#include "grins/physics_factory_variable_density_flow.h"

// GRINS
#include "grins/physics_factory_helper.h"
#include "grins/constant_viscosity.h"
#include "grins/constant_conductivity.h"
#include "grins/constant_specific_heat.h"

// Physics whose factories we're instantiating
#include "grins/low_mach_navier_stokes.h"
#include "grins/low_mach_navier_stokes_braack_stab.h"
#include "grins/low_mach_navier_stokes_spgsm_stab.h"
#include "grins/low_mach_navier_stokes_vms_stab.h"

namespace GRINS
{
  template<template<typename,typename,typename> class DerivedPhysics>
  libMesh::UniquePtr<Physics> PhysicsFactoryVariableDensityFlow<DerivedPhysics>::build_physics( const GetPot& input,
                                                                                                const std::string& physics_name )
  {
    std::string core_physics = this->find_core_physics_name(physics_name);

    std::string conductivity;
    PhysicsFactoryHelper::parse_conductivity_model(input,core_physics,conductivity);

    std::string viscosity;
    PhysicsFactoryHelper::parse_viscosity_model(input,core_physics,viscosity);

    std::string specific_heat;
    PhysicsFactoryHelper::parse_specific_heat_model(input,core_physics,specific_heat);

    libMesh::UniquePtr<Physics> new_physics;

    if(  conductivity == "constant" && viscosity == "constant" && specific_heat == "constant" )
      new_physics.reset( new DerivedPhysics<ConstantViscosity,ConstantSpecificHeat,ConstantConductivity>
                         (physics_name,input) );

    else
      this->prop_error_msg(physics_name, conductivity, viscosity, specific_heat);

    libmesh_assert(new_physics);

    return new_physics;
  }

  template<template<typename,typename,typename> class DerivedPhysics>
  void PhysicsFactoryVariableDensityFlow<DerivedPhysics>::prop_error_msg( const std::string& physics,
                                                                          const std::string& conductivity,
                                                                          const std::string& viscosity,
                                                                          const std::string& specific_heat ) const
  {
    std::string error = "================================================================\n";
    error += "Invalid combination of models for "+physics+"\n";
    error += "Viscosity model     = "+viscosity+"\n";
    error += "Conductivity model  = "+conductivity+"\n";
    error += "Specific heat model = "+specific_heat+"\n";
    error += "================================================================\n";

    libmesh_error_msg(error);
  }

  // Instantiate all the "variable density flow" Physics factories.
  PhysicsFactoryVariableDensityFlow<LowMachNavierStokes> grins_factory_low_mach_navier_stokes
  (PhysicsNaming::low_mach_navier_stokes(),PhysicsNaming::low_mach_navier_stokes());

  PhysicsFactoryVariableDensityFlow<LowMachNavierStokesSPGSMStabilization> grins_factory_lmns_spgsm_stab
  (PhysicsNaming::low_mach_navier_stokes_spgsm_stab(),PhysicsNaming::low_mach_navier_stokes());

  PhysicsFactoryVariableDensityFlow<LowMachNavierStokesVMSStabilization> grins_factory_lmns_vms_stab
  (PhysicsNaming::low_mach_navier_stokes_vms_stab(),PhysicsNaming::low_mach_navier_stokes());

  PhysicsFactoryVariableDensityFlow<LowMachNavierStokesBraackStabilization> grins_factory_lmns_braack_stab
  (PhysicsNaming::low_mach_navier_stokes_braack_stab(),PhysicsNaming::low_mach_navier_stokes());

} // end namespace GRINS
