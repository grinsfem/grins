//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_PHYSICS_NAMES_H
#define GRINS_PHYSICS_NAMES_H

#include <string>

namespace GRINS
{
  typedef std::string PhysicsName;

  const PhysicsName stokes = "Stokes";
  const PhysicsName incompressible_navier_stokes = "IncompressibleNavierStokes";
  const PhysicsName incompressible_navier_stokes_adjoint_stab = 
    "IncompressibleNavierStokesAdjointStabilization";
  const PhysicsName incompressible_navier_stokes_spgsm_stab = 
    "IncompressibleNavierStokesSPGSMStabilization";
  const PhysicsName velocity_drag = "VelocityDrag";
  const PhysicsName velocity_drag_adjoint_stab = "VelocityDragAdjointStabilization";
  const PhysicsName velocity_penalty = "VelocityPenalty";
  const PhysicsName velocity_penalty_adjoint_stab = "VelocityPenaltyAdjointStabilization";
  const PhysicsName averaged_fan = "AveragedFan";
  const PhysicsName averaged_turbine = "AveragedTurbine";
  const PhysicsName scalar_ode = "ScalarODE";
  const PhysicsName heat_conduction = "HeatConduction";
  const PhysicsName heat_transfer = "HeatTransfer";
  const PhysicsName heat_transfer_source = "HeatTransferSource";
  const PhysicsName heat_transfer_adjoint_stab = "HeatTransferAdjointStabilization";
  const PhysicsName heat_transfer_spgsm_stab = "HeatTransferSPGSMStabilization";
  const PhysicsName axisymmetric_heat_transfer = "AxisymmetricHeatTransfer";
  const PhysicsName boussinesq_buoyancy = "BoussinesqBuoyancy";
  const PhysicsName boussinesq_buoyancy_adjoint_stab = "BoussinesqBuoyancyAdjointStabilization";
  const PhysicsName boussinesq_buoyancy_spgsm_stab = "BoussinesqBuoyancySPGSMStabilization";
  const PhysicsName axisymmetric_boussinesq_buoyancy = "AxisymmetricBoussinesqBuoyancy";
  const PhysicsName low_mach_navier_stokes = "LowMachNavierStokes";
  const PhysicsName low_mach_navier_stokes_braack_stab = "LowMachNavierStokesBraackStabilization";
  const PhysicsName low_mach_navier_stokes_spgsm_stab = "LowMachNavierStokesSPGSMStabilization";
  const PhysicsName low_mach_navier_stokes_vms_stab = "LowMachNavierStokesVMSStabilization";
  const PhysicsName reacting_low_mach_navier_stokes = "ReactingLowMachNavierStokes";
  const PhysicsName elastic_membrane = "ElasticMembrane";
  const PhysicsName elastic_membrane_constant_pressure = "ElasticMembraneConstantPressure";
}

#endif //GRINS_PHYSICS_NAMES_H
