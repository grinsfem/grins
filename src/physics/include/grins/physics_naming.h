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


#ifndef GRINS_PHYSICS_NAMING_H
#define GRINS_PHYSICS_NAMING_H

#include <string>

namespace GRINS
{
  typedef std::string PhysicsName;

  class PhysicsNaming
  {
  public:

    static void set_suffix( const std::string& suff )
    { PhysicsNaming::suffix() = suff; }

    static void clear_suffix()
    { PhysicsNaming::suffix().clear(); }

    //! Extract the physics name and the suffix from the full_name
    /*! Note returned suffix includes the delimiter */
    static void extract_physics_and_suffix( const std::string& full_name,
                                            std::string& physics_name,
                                            std::string& suffix );

    //! Extract the physics name from the full_name.
    /*! If the delimiter is not present in the full_name, then
      the physics_name is the full_name. */
    static std::string extract_physics( const std::string& full_name );

    //! Extract the suffix from the full_name.
    /*! Note returned suffix includes the delimiter. If there is
      no delimiter present, then the returned suffix is an empty
      string. */
    static std::string extract_suffix( const std::string& full_name );

    static PhysicsName stokes()
    { return "Stokes"+suffix(); }

    static PhysicsName incompressible_navier_stokes()
    { return "IncompressibleNavierStokes"+suffix(); }

    static PhysicsName incompressible_navier_stokes_adjoint_stab()
    { return  "IncompressibleNavierStokesAdjointStabilization"+suffix(); }

    static PhysicsName incompressible_navier_stokes_spgsm_stab()
    { return "IncompressibleNavierStokesSPGSMStabilization"+suffix(); }

    static PhysicsName velocity_drag()
    { return "VelocityDrag"+suffix(); }

    static PhysicsName velocity_drag_adjoint_stab()
    { return "VelocityDragAdjointStabilization"+suffix(); }

    static PhysicsName velocity_penalty()
    { return "VelocityPenalty"+suffix(); }

    static PhysicsName velocity_penalty2()
    { return "VelocityPenalty2"+suffix(); }

    static PhysicsName velocity_penalty3()
    { return "VelocityPenalty3"+suffix(); }

    static PhysicsName velocity_penalty_adjoint_stab()
    { return "VelocityPenaltyAdjointStabilization"+suffix(); }

    static PhysicsName velocity_penalty2_adjoint_stab()
    { return "VelocityPenalty2AdjointStabilization"+suffix(); }

    static PhysicsName velocity_penalty3_adjoint_stab()
    { return "VelocityPenalty3AdjointStabilization"+suffix(); }

    static PhysicsName parsed_velocity_source()
    { return "ParsedVelocitySource"+suffix(); }

    static PhysicsName parsed_velocity_source_adjoint_stab()
    { return "ParsedVelocitySourceAdjointStabilization"+suffix(); }

    static PhysicsName averaged_fan()
    { return "AveragedFan"+suffix(); }

    static PhysicsName averaged_fan_adjoint_stab()
    { return "AveragedFanAdjointStabilization"+suffix(); }

    static PhysicsName averaged_turbine()
    { return "AveragedTurbine"+suffix(); }

    static PhysicsName spalart_allmaras()
    { return "SpalartAllmaras"+suffix(); }

    static PhysicsName spalart_allmaras_spgsm_stab()
    { return "SpalartAllmarasSPGSMStabilization"+suffix(); }

    static PhysicsName scalar_ode()
    { return "ScalarODE"+suffix(); }

    static PhysicsName heat_conduction()
    { return "HeatConduction"+suffix(); }

    static PhysicsName heat_transfer()
    { return "HeatTransfer"+suffix(); }

    static PhysicsName heat_transfer_adjoint_stab()
    { return "HeatTransferAdjointStabilization"+suffix(); }

    static PhysicsName heat_transfer_spgsm_stab()
    { return "HeatTransferSPGSMStabilization"+suffix(); }

    static PhysicsName axisymmetric_heat_transfer()
    { return "AxisymmetricHeatTransfer"+suffix(); }

    static PhysicsName boussinesq_buoyancy()
    { return "BoussinesqBuoyancy"+suffix(); }

    static PhysicsName boussinesq_buoyancy_adjoint_stab()
    { return "BoussinesqBuoyancyAdjointStabilization"+suffix(); }

    static PhysicsName boussinesq_buoyancy_spgsm_stab()
    { return "BoussinesqBuoyancySPGSMStabilization"+suffix(); }

    static PhysicsName axisymmetric_boussinesq_buoyancy()
    { return "AxisymmetricBoussinesqBuoyancy"+suffix(); }

    static PhysicsName low_mach_navier_stokes()
    { return "LowMachNavierStokes"+suffix(); }

    static PhysicsName low_mach_navier_stokes_braack_stab()
    { return "LowMachNavierStokesBraackStabilization"+suffix(); }

    static PhysicsName low_mach_navier_stokes_spgsm_stab()
    { return "LowMachNavierStokesSPGSMStabilization"+suffix(); }

    static PhysicsName low_mach_navier_stokes_vms_stab()
    { return "LowMachNavierStokesVMSStabilization"+suffix(); }

    static PhysicsName reacting_low_mach_navier_stokes()
    { return "ReactingLowMachNavierStokes"+suffix(); }

    static PhysicsName reacting_low_mach_navier_stokes_spgsm_stab()
    { return "ReactingLowMachNavierStokesSPGSMStabilization"+suffix(); }

    static PhysicsName od_premixed_flame()
    { return "ODPremixedFlame"+suffix(); }

    static PhysicsName elastic_membrane()
    { return "ElasticMembrane"+suffix(); }

    static PhysicsName elastic_cable()
    { return "ElasticCable"+suffix(); }

    static PhysicsName elastic_cable_rayleigh_damping()
    { return "ElasticCableRayleighDamping"+suffix(); }

    static PhysicsName elastic_membrane_rayleigh_damping()
    { return "ElasticMembraneRayleighDamping"+suffix(); }

    static PhysicsName elastic_membrane_constant_pressure()
    { return "ElasticMembraneConstantPressure"+suffix(); }

    static PhysicsName elastic_cable_constant_gravity()
    { return "ElasticCableConstantGravity"+suffix(); }

    static PhysicsName constant_source_term()
    { return "ConstantSourceTerm"+suffix(); }

    static PhysicsName parsed_source_term()
    { return "ParsedSourceTerm"+suffix(); }

    static PhysicsName convection_diffusion()
    { return "ConvectionDiffusion"+suffix(); }

    static PhysicsName variable_pinning()
    { return "VariablePinning"+suffix(); }

  private:

    static std::string physics_name_delimiter()
    { return ":"; }

    static std::string& suffix()
    {
      static std::string _suffix;
      return _suffix;
    };

  };

}

#endif //GRINS_PHYSICS_NAMING_H
