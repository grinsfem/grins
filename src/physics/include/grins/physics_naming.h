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


#ifndef GRINS_PHYSICS_NAMING_H
#define GRINS_PHYSICS_NAMING_H

#include <string>

namespace GRINS
{
  typedef std::string PhysicsName;

  class PhysicsNaming
  {
  public:

    static PhysicsName stokes()
    { return "Stokes"; }

    static PhysicsName incompressible_navier_stokes()
    { return "IncompressibleNavierStokes"; }

    static PhysicsName incompressible_navier_stokes_adjoint_stab()
    { return  "IncompressibleNavierStokesAdjointStabilization"; }

    static PhysicsName incompressible_navier_stokes_spgsm_stab()
    { return "IncompressibleNavierStokesSPGSMStabilization"; }

    static PhysicsName velocity_drag()
    { return "VelocityDrag"; }

    static PhysicsName velocity_drag_adjoint_stab()
    { return "VelocityDragAdjointStabilization"; }

    static PhysicsName velocity_penalty()
    { return "VelocityPenalty"; }

    static PhysicsName velocity_penalty2()
    { return "VelocityPenalty2"; }

    static PhysicsName velocity_penalty3()
    { return "VelocityPenalty3"; }

    static PhysicsName velocity_penalty_adjoint_stab()
    { return "VelocityPenaltyAdjointStabilization"; }

    static PhysicsName velocity_penalty2_adjoint_stab()
    { return "VelocityPenalty2AdjointStabilization"; }

    static PhysicsName velocity_penalty3_adjoint_stab()
    { return "VelocityPenalty3AdjointStabilization"; }

    static PhysicsName parsed_velocity_source()
    { return "ParsedVelocitySource"; }

    static PhysicsName parsed_velocity_source_adjoint_stab()
    { return "ParsedVelocitySourceAdjointStabilization"; }

    static PhysicsName averaged_fan()
    { return "AveragedFan"; }

    static PhysicsName averaged_fan_adjoint_stab()
    { return "AveragedFanAdjointStabilization"; }

    static PhysicsName averaged_turbine()
    { return "AveragedTurbine"; }

    static PhysicsName spalart_allmaras()
    { return "SpalartAllmaras"; }

    static PhysicsName spalart_allmaras_spgsm_stab()
    { return "SpalartAllmarasSPGSMStabilization"; }

    static PhysicsName scalar_ode()
    { return "ScalarODE"; }

    static PhysicsName heat_conduction()
    { return "HeatConduction"; }

    static PhysicsName heat_transfer()
    { return "HeatTransfer"; }

    static PhysicsName heat_transfer_source()
    { return "HeatTransferSource"; }

    static PhysicsName heat_transfer_adjoint_stab()
    { return "HeatTransferAdjointStabilization"; }

    static PhysicsName heat_transfer_spgsm_stab()
    { return "HeatTransferSPGSMStabilization"; }

    static PhysicsName axisymmetric_heat_transfer()
    { return "AxisymmetricHeatTransfer"; }

    static PhysicsName boussinesq_buoyancy()
    { return "BoussinesqBuoyancy"; }

    static PhysicsName boussinesq_buoyancy_adjoint_stab()
    { return "BoussinesqBuoyancyAdjointStabilization"; }

    static PhysicsName boussinesq_buoyancy_spgsm_stab()
    { return "BoussinesqBuoyancySPGSMStabilization"; }

    static PhysicsName axisymmetric_boussinesq_buoyancy()
    { return "AxisymmetricBoussinesqBuoyancy"; }

    static PhysicsName low_mach_navier_stokes()
    { return "LowMachNavierStokes"; }

    static PhysicsName low_mach_navier_stokes_braack_stab()
    { return "LowMachNavierStokesBraackStabilization"; }

    static PhysicsName low_mach_navier_stokes_spgsm_stab()
    { return "LowMachNavierStokesSPGSMStabilization"; }

    static PhysicsName low_mach_navier_stokes_vms_stab()
    { return "LowMachNavierStokesVMSStabilization"; }

    static PhysicsName reacting_low_mach_navier_stokes()
    { return "ReactingLowMachNavierStokes"; }

    static PhysicsName elastic_membrane()
    { return "ElasticMembrane"; }

    static PhysicsName elastic_cable()
    { return "ElasticCable"; }

    static PhysicsName elastic_membrane_constant_pressure()
    { return "ElasticMembraneConstantPressure"; }

    static PhysicsName elastic_cable_constant_gravity()
    { return "ElasticCableConstantGravity"; }

    static PhysicsName constant_source_term()
    { return "ConstantSourceTerm"; }

    static PhysicsName parsed_source_term()
    { return "ParsedSourceTerm"; }
  };

}

#endif //GRINS_PHYSICS_NAMING_H
