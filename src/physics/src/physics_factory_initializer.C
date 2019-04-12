//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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
#include "grins/physics_factory_initializer.h"

// GRINS
#include "grins/physics_factory_basic.h"
#include "grins/physics_factory_heat_transfer.h"
#include "grins/physics_factory_incompressible_flow.h"
#include "grins/physics_factory_incompressible_turb_flow.h"
#include "grins/physics_factory_one_d_stress_solids.h"
#include "grins/physics_factory_plane_stress_solids.h"
#include "grins/physics_factory_variable_density_flow.h"
#include "grins/physics_factory_reacting_flows.h"

#include "grins/physics_naming.h"
#include "grins/scalar_ode.h"
#include "grins/boussinesq_buoyancy.h"
#include "grins/axisym_boussinesq_buoyancy.h"
#include "grins/constant_source_term.h"
#include "grins/parsed_source_term.h"
#include "grins/convection_diffusion.h"
#include "grins/variable_pinning.h"

#include "grins/heat_conduction.h"
#include "grins/heat_transfer.h"
#include "grins/heat_transfer_adjoint_stab.h"
#include "grins/heat_transfer_spgsm_stab.h"
#include "grins/axisym_heat_transfer.h"

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

#include "grins/spalart_allmaras.h"
#include "grins/spalart_allmaras_spgsm_stab.h"

#include "grins/elastic_cable.h"
#include "grins/elastic_cable_rayleigh_damping.h"

#include "grins/elastic_membrane.h"
#include "grins/elastic_membrane_rayleigh_damping.h"
#include "grins/elastic_membrane_constant_pressure.h"
#include "grins/elastic_membrane_pressure.h"
#include "grins/parsed_pressure.h"

#include "grins/low_mach_navier_stokes.h"
#include "grins/low_mach_navier_stokes_braack_stab.h"
#include "grins/low_mach_navier_stokes_spgsm_stab.h"
#include "grins/low_mach_navier_stokes_vms_stab.h"

namespace GRINS
{
  PhysicsFactoryInitializer::PhysicsFactoryInitializer()
  {
    static PhysicsFactoryBasic<ScalarODE> grins_factory_scalar_ode(PhysicsNaming::scalar_ode());

    static PhysicsFactoryBasic<BoussinesqBuoyancy> grins_factory_boussinesq(PhysicsNaming::boussinesq_buoyancy());

    // This one needs to die. Regular Boussinesq should handle the axisymmetry
    static PhysicsFactoryBasic<AxisymmetricBoussinesqBuoyancy>
      grins_factory_axi_boussinesq(PhysicsNaming::axisymmetric_boussinesq_buoyancy());

    static PhysicsFactoryBasic<ElasticMembraneConstantPressure>
      grins_factory_elastic_membrane_constant_pressure(PhysicsNaming::elastic_membrane_constant_pressure());

    static PhysicsFactoryBasic<ElasticMembranePressure<ParsedPressure> >
      grins_factory_elastic_membrane_parsed_pressure(PhysicsNaming::elastic_membrane_parsed_pressure());

    // This one needs to die and just have the parsed version
    static PhysicsFactoryBasic<ConstantSourceTerm>
      grins_factory_constant_source_term(PhysicsNaming::constant_source_term());

    static PhysicsFactoryBasic<ParsedSourceTerm>
      grins_factory_parsed_source_term(PhysicsNaming::parsed_source_term());

    static PhysicsFactoryBasic<ConvectionDiffusion>
      grins_factory_convection_difffusion(PhysicsNaming::convection_diffusion());

    static PhysicsFactoryBasic<VariablePinning>
      grins_factory_variable_pinning(PhysicsNaming::variable_pinning());


    static PhysicsFactoryHeatTransfer<HeatConduction> grins_factory_heat_conduction
      (PhysicsNaming::heat_conduction(),PhysicsNaming::heat_conduction());

    static PhysicsFactoryHeatTransfer<HeatTransfer> grins_factory_heat_transfer
      (PhysicsNaming::heat_transfer(),PhysicsNaming::heat_transfer());

    static PhysicsFactoryHeatTransfer<HeatTransferAdjointStabilization> grins_factory_heat_transfer_adjoint_stab
      (PhysicsNaming::heat_transfer_adjoint_stab(),PhysicsNaming::heat_transfer());

    static PhysicsFactoryHeatTransfer<HeatTransferAdjointStabilization> grins_factory_heat_transfer_spgsm_stab
      (PhysicsNaming::heat_transfer_spgsm_stab(),PhysicsNaming::heat_transfer());

    // This needs to die. Axisymmetry should be handled within heat_transfer
    static PhysicsFactoryHeatTransfer<AxisymmetricHeatTransfer> grins_factory_axi_heat_transfer
      (PhysicsNaming::axisymmetric_heat_transfer(),PhysicsNaming::axisymmetric_heat_transfer());


    static PhysicsFactoryIncompressibleFlow<IncompressibleNavierStokes> grins_factory_incompressible_navier_stokes
      (PhysicsNaming::incompressible_navier_stokes(),PhysicsNaming::incompressible_navier_stokes());

    static PhysicsFactoryIncompressibleFlow<Stokes> grins_factory_stokes
      (PhysicsNaming::stokes(),PhysicsNaming::stokes());

    static PhysicsFactoryIncompressibleFlow<IncompressibleNavierStokesAdjointStabilization>
      grins_factory_ins_adjoint_stab
      (PhysicsNaming::incompressible_navier_stokes_adjoint_stab(),PhysicsNaming::incompressible_navier_stokes());

    static PhysicsFactoryIncompressibleFlow<IncompressibleNavierStokesSPGSMStabilization>
      grins_factory_spsgm_adjoint_stab
      (PhysicsNaming::incompressible_navier_stokes_spgsm_stab(),PhysicsNaming::incompressible_navier_stokes());

    static PhysicsFactoryIncompressibleFlow<VelocityDrag> grins_factory_velocity_drag
      (PhysicsNaming::velocity_drag(),PhysicsNaming::incompressible_navier_stokes());

    static PhysicsFactoryIncompressibleFlow<VelocityDragAdjointStabilization>
      grins_factory_velocity_drag_adjoint_stab
      (PhysicsNaming::velocity_drag_adjoint_stab(),PhysicsNaming::incompressible_navier_stokes());

    static PhysicsFactoryIncompressibleFlow<VelocityPenalty> grins_factory_vel_penalty
      (PhysicsNaming::velocity_penalty(),PhysicsNaming::incompressible_navier_stokes());

    static PhysicsFactoryIncompressibleFlow<VelocityPenalty> grins_factory_vel_penalty2
      (PhysicsNaming::velocity_penalty2(),PhysicsNaming::incompressible_navier_stokes());

    static PhysicsFactoryIncompressibleFlow<VelocityPenalty> grins_factory_vel_penalty3
      (PhysicsNaming::velocity_penalty3(),PhysicsNaming::incompressible_navier_stokes());

    static PhysicsFactoryIncompressibleFlow<VelocityPenaltyAdjointStabilization>
      grins_factory_vel_penalty_adjoint_stab
      (PhysicsNaming::velocity_penalty_adjoint_stab(),PhysicsNaming::incompressible_navier_stokes());

    static PhysicsFactoryIncompressibleFlow<VelocityPenaltyAdjointStabilization>
      grins_factory_vel_penalty2_adjoint_stab
      (PhysicsNaming::velocity_penalty2_adjoint_stab(),PhysicsNaming::incompressible_navier_stokes());

    static PhysicsFactoryIncompressibleFlow<VelocityPenaltyAdjointStabilization>
      grins_factory_vel_penalty3_adjoint_stab
      (PhysicsNaming::velocity_penalty3_adjoint_stab(),PhysicsNaming::incompressible_navier_stokes());

    static PhysicsFactoryIncompressibleFlow<ParsedVelocitySource> grins_factory_parsed_vel_source
      (PhysicsNaming::parsed_velocity_source(),PhysicsNaming::incompressible_navier_stokes());

    static PhysicsFactoryIncompressibleFlow<ParsedVelocitySourceAdjointStabilization>
      grins_factory_parsed_vel_source_adjoint_stab
      (PhysicsNaming::parsed_velocity_source_adjoint_stab(),PhysicsNaming::incompressible_navier_stokes());

    static PhysicsFactoryIncompressibleFlow<AveragedFan> grins_factory_averaged_fan
      (PhysicsNaming::averaged_fan(),PhysicsNaming::incompressible_navier_stokes());

    static PhysicsFactoryIncompressibleFlow<AveragedFanAdjointStabilization>
      grins_factory_averaged_fan_adjoint_stab
      (PhysicsNaming::averaged_fan_adjoint_stab(),PhysicsNaming::incompressible_navier_stokes());

    static PhysicsFactoryIncompressibleFlow<AveragedTurbine> grins_factory_averaged_turbine
      (PhysicsNaming::averaged_turbine(),PhysicsNaming::incompressible_navier_stokes());

    static PhysicsFactoryIncompressibleFlow<BoussinesqBuoyancyAdjointStabilization>
      grins_factory_boussinesq_adjoint_stab
      (PhysicsNaming::boussinesq_buoyancy_adjoint_stab(),PhysicsNaming::boussinesq_buoyancy());

    static PhysicsFactoryIncompressibleFlow<BoussinesqBuoyancySPGSMStabilization>
      grins_factory_boussinesq_spgsm_stab
      (PhysicsNaming::boussinesq_buoyancy_spgsm_stab(),PhysicsNaming::boussinesq_buoyancy());


    static PhysicsFactoryIncompressibleTurbFlow<SpalartAllmaras> grins_factory_spalart_allmaras
      (PhysicsNaming::spalart_allmaras(),PhysicsNaming::spalart_allmaras());

    static PhysicsFactoryIncompressibleTurbFlow<SpalartAllmarasSPGSMStabilization>
      grins_factory_spalart_allmaras_spgsm_stab
      (PhysicsNaming::spalart_allmaras_spgsm_stab(),PhysicsNaming::spalart_allmaras());


    static PhysicsFactoryOneDStressSolids<ElasticCable> grins_factory_elastic_cable
      (PhysicsNaming::elastic_cable(),PhysicsNaming::elastic_cable());

    static PhysicsFactoryOneDStressSolids<ElasticCableRayleighDamping> grins_factory_elastic_cable_rayleigh_damping
      (PhysicsNaming::elastic_cable_rayleigh_damping(),PhysicsNaming::elastic_cable());


    static PhysicsFactoryPlaneStressSolids<ElasticMembrane> grins_factory_elastic_membrane
      (PhysicsNaming::elastic_membrane(),PhysicsNaming::elastic_membrane());

    static PhysicsFactoryPlaneStressSolids<ElasticMembraneRayleighDamping>
      grins_factory_elastic_membrane_rayleigh_damping
      (PhysicsNaming::elastic_membrane_rayleigh_damping(),PhysicsNaming::elastic_membrane());


    static PhysicsFactoryVariableDensityFlow<LowMachNavierStokes> grins_factory_low_mach_navier_stokes
      (PhysicsNaming::low_mach_navier_stokes(),PhysicsNaming::low_mach_navier_stokes());

    static PhysicsFactoryVariableDensityFlow<LowMachNavierStokesSPGSMStabilization> grins_factory_lmns_spgsm_stab
      (PhysicsNaming::low_mach_navier_stokes_spgsm_stab(),PhysicsNaming::low_mach_navier_stokes());

    static PhysicsFactoryVariableDensityFlow<LowMachNavierStokesVMSStabilization> grins_factory_lmns_vms_stab
      (PhysicsNaming::low_mach_navier_stokes_vms_stab(),PhysicsNaming::low_mach_navier_stokes());

    static PhysicsFactoryVariableDensityFlow<LowMachNavierStokesBraackStabilization> grins_factory_lmns_braack_stab
      (PhysicsNaming::low_mach_navier_stokes_braack_stab(),PhysicsNaming::low_mach_navier_stokes());

    ReactingFlowsPhysicsFactoryInitializer grins_factory_reacting_flows_init;
  }
}// end namespace GRINS
