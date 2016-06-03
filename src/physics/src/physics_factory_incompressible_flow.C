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
#include "grins/physics_factory_incompressible_flow.h"

// GRINS
#include "grins/physics_factory_helper.h"
#include "grins/constant_viscosity.h"
#include "grins/parsed_viscosity.h"
#include "grins/spalart_allmaras_viscosity.h"

// Physics we're instantiating
#include "grins/inc_navier_stokes.h"
#include "grins/stokes.h"
#include "grins/inc_navier_stokes_adjoint_stab.h"
#include "grins/inc_navier_stokes_spgsm_stab.h"
#include "grins/velocity_drag.h"
#include "grins/velocity_drag_adjoint_stab.h"
#include "grins/velocity_penalty.h"
#include "grins/velocity_penalty_adjoint_stab.h"
#include "grins/parsed_velocity_source.h"
#include "grins/parsed_velocity_source_adjoint_stab.h"
#include "grins/averaged_fan.h"
#include "grins/averaged_fan_adjoint_stab.h"
#include "grins/averaged_turbine.h"
#include "grins/boussinesq_buoyancy_adjoint_stab.h"
#include "grins/boussinesq_buoyancy_spgsm_stab.h"

namespace GRINS
{
  template<template<typename> class DerivedPhysics>
  libMesh::UniquePtr<Physics> PhysicsFactoryIncompressibleFlow<DerivedPhysics>::build_physics( const GetPot& input,
                                                                                               const std::string& physics_name )
  {
    std::string core_physics = this->find_core_physics_name(physics_name);

    std::string viscosity;
    PhysicsFactoryHelper::parse_viscosity_model(input,core_physics,viscosity);

    libMesh::UniquePtr<Physics> new_physics;

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
  void PhysicsFactoryIncompressibleFlow<DerivedPhysics>::visc_error_msg( const std::string& physics,
                                                                         const std::string& viscosity ) const
  {
    std::string error = "================================================================\n";
    error += "Invalid viscosity model for "+physics+"\n";
    error += "Viscosity model     = "+viscosity+"\n";
    error += "================================================================\n";

    libmesh_error_msg(error);
  }

  // Instantiate all the "incompressble flow" Physics factories.
  PhysicsFactoryIncompressibleFlow<IncompressibleNavierStokes> grins_factory_incompressible_navier_stokes
  (PhysicsNaming::incompressible_navier_stokes(),PhysicsNaming::incompressible_navier_stokes());

  PhysicsFactoryIncompressibleFlow<Stokes> grins_factory_stokes
  (PhysicsNaming::stokes(),PhysicsNaming::stokes());

  PhysicsFactoryIncompressibleFlow<IncompressibleNavierStokesAdjointStabilization> grins_factory_ins_adjoint_stab
  (PhysicsNaming::incompressible_navier_stokes_adjoint_stab(),PhysicsNaming::incompressible_navier_stokes());

  PhysicsFactoryIncompressibleFlow<IncompressibleNavierStokesSPGSMStabilization> grins_factory_spsgm_adjoint_stab
  (PhysicsNaming::incompressible_navier_stokes_spgsm_stab(),PhysicsNaming::incompressible_navier_stokes());

  PhysicsFactoryIncompressibleFlow<VelocityDrag> grins_factory_velocity_drag
  (PhysicsNaming::velocity_drag(),PhysicsNaming::incompressible_navier_stokes());

  PhysicsFactoryIncompressibleFlow<VelocityDragAdjointStabilization> grins_factory_velocity_drag_adjoint_stab
  (PhysicsNaming::velocity_drag_adjoint_stab(),PhysicsNaming::incompressible_navier_stokes());

  PhysicsFactoryIncompressibleFlow<VelocityPenalty> grins_factory_vel_penalty
  (PhysicsNaming::velocity_penalty(),PhysicsNaming::incompressible_navier_stokes());

  PhysicsFactoryIncompressibleFlow<VelocityPenalty> grins_factory_vel_penalty2
  (PhysicsNaming::velocity_penalty2(),PhysicsNaming::incompressible_navier_stokes());

  PhysicsFactoryIncompressibleFlow<VelocityPenalty> grins_factory_vel_penalty3
  (PhysicsNaming::velocity_penalty3(),PhysicsNaming::incompressible_navier_stokes());

  PhysicsFactoryIncompressibleFlow<VelocityPenaltyAdjointStabilization> grins_factory_vel_penalty_adjoint_stab
  (PhysicsNaming::velocity_penalty_adjoint_stab(),PhysicsNaming::incompressible_navier_stokes());

  PhysicsFactoryIncompressibleFlow<VelocityPenaltyAdjointStabilization> grins_factory_vel_penalty2_adjoint_stab
  (PhysicsNaming::velocity_penalty2_adjoint_stab(),PhysicsNaming::incompressible_navier_stokes());

  PhysicsFactoryIncompressibleFlow<VelocityPenaltyAdjointStabilization> grins_factory_vel_penalty3_adjoint_stab
  (PhysicsNaming::velocity_penalty3_adjoint_stab(),PhysicsNaming::incompressible_navier_stokes());

  PhysicsFactoryIncompressibleFlow<ParsedVelocitySource> grins_factory_parsed_vel_source
  (PhysicsNaming::parsed_velocity_source(),PhysicsNaming::incompressible_navier_stokes());

  PhysicsFactoryIncompressibleFlow<ParsedVelocitySourceAdjointStabilization> grins_factory_parsed_vel_source_adjoint_stab
  (PhysicsNaming::parsed_velocity_source_adjoint_stab(),PhysicsNaming::incompressible_navier_stokes());

  PhysicsFactoryIncompressibleFlow<AveragedFan> grins_factory_averaged_fan
  (PhysicsNaming::averaged_fan(),PhysicsNaming::incompressible_navier_stokes());

  PhysicsFactoryIncompressibleFlow<AveragedFanAdjointStabilization> grins_factory_averaged_fan_adjoint_stab
  (PhysicsNaming::averaged_fan_adjoint_stab(),PhysicsNaming::incompressible_navier_stokes());

  PhysicsFactoryIncompressibleFlow<AveragedTurbine> grins_factory_averaged_turbine
  (PhysicsNaming::averaged_turbine(),PhysicsNaming::incompressible_navier_stokes());

  PhysicsFactoryIncompressibleFlow<BoussinesqBuoyancyAdjointStabilization> grins_factory_boussinesq_adjoint_stab
  (PhysicsNaming::boussinesq_buoyancy_adjoint_stab(),PhysicsNaming::boussinesq_buoyancy());

  PhysicsFactoryIncompressibleFlow<BoussinesqBuoyancySPGSMStabilization> grins_factory_boussinesq_spgsm_stab
  (PhysicsNaming::boussinesq_buoyancy_spgsm_stab(),PhysicsNaming::boussinesq_buoyancy());

} // end namespace GRINS
