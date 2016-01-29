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

    static void set_suffix( const std::string& suffix )
    { _suffix = suffix; }

    static void clear_suffix()
    { _suffix.clear(); }

    static PhysicsName stokes()
    { return "Stokes"+_suffix; }

    static PhysicsName incompressible_navier_stokes()
    { return "IncompressibleNavierStokes"+_suffix; }

    static PhysicsName incompressible_navier_stokes_adjoint_stab()
    { return  "IncompressibleNavierStokesAdjointStabilization"+_suffix; }

    static PhysicsName incompressible_navier_stokes_spgsm_stab()
    { return "IncompressibleNavierStokesSPGSMStabilization"+_suffix; }

    static PhysicsName velocity_drag()
    { return "VelocityDrag"+_suffix; }

    static PhysicsName velocity_drag_adjoint_stab()
    { return "VelocityDragAdjointStabilization"+_suffix; }

    static PhysicsName velocity_penalty()
    { return "VelocityPenalty"+_suffix; }

    static PhysicsName velocity_penalty2()
    { return "VelocityPenalty2"+_suffix; }

    static PhysicsName velocity_penalty3()
    { return "VelocityPenalty3"+_suffix; }

    static PhysicsName velocity_penalty_adjoint_stab()
    { return "VelocityPenaltyAdjointStabilization"+_suffix; }

    static PhysicsName velocity_penalty2_adjoint_stab()
    { return "VelocityPenalty2AdjointStabilization"+_suffix; }

    static PhysicsName velocity_penalty3_adjoint_stab()
    { return "VelocityPenalty3AdjointStabilization"+_suffix; }

    static PhysicsName parsed_velocity_source()
    { return "ParsedVelocitySource"+_suffix; }

    static PhysicsName parsed_velocity_source_adjoint_stab()
    { return "ParsedVelocitySourceAdjointStabilization"+_suffix; }

    static PhysicsName averaged_fan()
    { return "AveragedFan"+_suffix; }

    static PhysicsName averaged_fan_adjoint_stab()
    { return "AveragedFanAdjointStabilization"+_suffix; }

    static PhysicsName averaged_turbine()
    { return "AveragedTurbine"+_suffix; }

    static PhysicsName spalart_allmaras()
    { return "SpalartAllmaras"+_suffix; }

    static PhysicsName spalart_allmaras_spgsm_stab()
    { return "SpalartAllmarasSPGSMStabilization"+_suffix; }

    static PhysicsName scalar_ode()
    { return "ScalarODE"+_suffix; }

    static PhysicsName heat_conduction()
    { return "HeatConduction"+_suffix; }

    static PhysicsName heat_transfer()
    { return "HeatTransfer"+_suffix; }

    static PhysicsName heat_transfer_source()
    { return "HeatTransferSource"+_suffix; }

    static PhysicsName heat_transfer_adjoint_stab()
    { return "HeatTransferAdjointStabilization"+_suffix; }

    static PhysicsName heat_transfer_spgsm_stab()
    { return "HeatTransferSPGSMStabilization"+_suffix; }

    static PhysicsName axisymmetric_heat_transfer()
    { return "AxisymmetricHeatTransfer"+_suffix; }

    static PhysicsName boussinesq_buoyancy()
    { return "BoussinesqBuoyancy"+_suffix; }

    static PhysicsName boussinesq_buoyancy_adjoint_stab()
    { return "BoussinesqBuoyancyAdjointStabilization"+_suffix; }

    static PhysicsName boussinesq_buoyancy_spgsm_stab()
    { return "BoussinesqBuoyancySPGSMStabilization"+_suffix; }

    static PhysicsName axisymmetric_boussinesq_buoyancy()
    { return "AxisymmetricBoussinesqBuoyancy"+_suffix; }

    static PhysicsName low_mach_navier_stokes()
    { return "LowMachNavierStokes"+_suffix; }

    static PhysicsName low_mach_navier_stokes_braack_stab()
    { return "LowMachNavierStokesBraackStabilization"+_suffix; }

    static PhysicsName low_mach_navier_stokes_spgsm_stab()
    { return "LowMachNavierStokesSPGSMStabilization"+_suffix; }

    static PhysicsName low_mach_navier_stokes_vms_stab()
    { return "LowMachNavierStokesVMSStabilization"+_suffix; }

    static PhysicsName reacting_low_mach_navier_stokes()
    { return "ReactingLowMachNavierStokes"+_suffix; }

    static PhysicsName elastic_membrane()
    { return "ElasticMembrane"+_suffix; }

    static PhysicsName elastic_cable()
    { return "ElasticCable"+_suffix; }

    static PhysicsName elastic_membrane_constant_pressure()
    { return "ElasticMembraneConstantPressure"+_suffix; }

    static PhysicsName elastic_cable_constant_gravity()
    { return "ElasticCableConstantGravity"+_suffix; }

    static PhysicsName constant_source_term()
    { return "ConstantSourceTerm"+_suffix; }

    static PhysicsName parsed_source_term()
    { return "ParsedSourceTerm"+_suffix; }

  private:

    static std::string _suffix;

  };

}

#endif //GRINS_PHYSICS_NAMING_H
