// This class
#include "grins/physics_factory_reacting_flows.h"

// GRINS
#include "grins/physics_factory_helper.h"
#include "grins/cantera_mixture.h"
#include "grins/cantera_evaluator.h"
#include "grins/antioch_mixture_averaged_transport_mixture.h"
#include "grins/antioch_mixture_averaged_transport_evaluator.h"
#include "grins/antioch_constant_transport_mixture.h"
#include "grins/antioch_constant_transport_evaluator.h"

// Antioch
#ifdef GRINS_HAVE_ANTIOCH
#include "antioch_config.h"
#include "antioch/sutherland_viscosity.h"
#include "antioch/blottner_viscosity.h"
#include "antioch/sutherland_parsing.h"
#include "antioch/blottner_parsing.h"
#include "antioch/eucken_thermal_conductivity.h"
#include "antioch/constant_lewis_diffusivity.h"
#endif // GRINS_HAVE_ANTIOCH

// Physics whose factories we're instantiating
#include "grins/reacting_low_mach_navier_stokes.h"
#include "grins/reacting_low_mach_navier_stokes_spgsm_stab.h"

namespace GRINS
{
  template<template<typename,typename> class DerivedPhysics>
  libMesh::UniquePtr<Physics> PhysicsFactoryReactingFlows<DerivedPhysics>::build_physics( const GetPot& input,
                                                                                          const std::string& physics_name )
  {
    std::string core_physics = this->find_core_physics_name(physics_name);

    std::string thermochem_lib;
    PhysicsFactoryHelper::parse_thermochemistry_model( input,
                                                       core_physics,
                                                       thermochem_lib );

    libMesh::UniquePtr<Physics> new_physics;

    if( thermochem_lib == "cantera" )
      {
#ifdef GRINS_HAVE_CANTERA
        new_physics.reset(new DerivedPhysics<CanteraMixture,CanteraEvaluator>(physics_name,input));
#else
        libmesh_error_msg("Error: Cantera not enabled in this configuration. Reconfigure using --with-cantera option.");

#endif // GRINS_HAVE_CANTERA
      }

    else if( thermochem_lib == "antioch" )
      {
#ifdef GRINS_HAVE_ANTIOCH

        std::string transport_model;
        std::string thermo_model;
        std::string viscosity_model;
        std::string conductivity_model;
        std::string diffusivity_model;

        PhysicsFactoryHelper::parse_antioch_models( input,
                                                    core_physics,
                                                    transport_model,
                                                    thermo_model,
                                                    viscosity_model,
                                                    conductivity_model,
                                                    diffusivity_model );

        if( transport_model == std::string("mixture_averaged") )
          {
            if( (thermo_model == std::string("stat_mech")) &&
                (diffusivity_model == std::string("constant_lewis")) &&
                (conductivity_model == std::string("eucken")) &&
                (viscosity_model == std::string("sutherland")) )
              {
                new_physics.reset(new DerivedPhysics<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                                   Antioch::SutherlandViscosity<libMesh::Real>,
                                                                                                   Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                                                   Antioch::ConstantLewisDiffusivity<libMesh::Real> >,
                                                     GRINS::AntiochMixtureAveragedTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                                     Antioch::SutherlandViscosity<libMesh::Real>,
                                                                                                     Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                                                     Antioch::ConstantLewisDiffusivity<libMesh::Real> > >
                                  (physics_name,input) );
              }
            else if( (thermo_model == std::string("stat_mech")) &&
                     (diffusivity_model == std::string("constant_lewis")) &&
                     (conductivity_model == std::string("eucken")) &&
                     (viscosity_model == std::string("blottner")) )
              {
                new_physics.reset(new DerivedPhysics<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                                   Antioch::BlottnerViscosity<libMesh::Real>,
                                                                                                   Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                                                   Antioch::ConstantLewisDiffusivity<libMesh::Real> >,
                                                     GRINS::AntiochMixtureAveragedTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                                     Antioch::BlottnerViscosity<libMesh::Real>,
                                                                                                     Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                                                     Antioch::ConstantLewisDiffusivity<libMesh::Real> > >(physics_name,input) );
              }
            else if( (thermo_model == std::string("stat_mech")) &&
                     (diffusivity_model == std::string("kinetics_theory")) &&
                     (conductivity_model == std::string("kinetics_theory")) &&
                     (viscosity_model == std::string("kinetics_theory")) )
              {
#ifdef ANTIOCH_HAVE_GSL
                new_physics.reset(new DerivedPhysics<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                                   Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>,
                                                                                                   Antioch::KineticsTheoryThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real>,libMesh::Real>,
                                                                                                   Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> >,
                                                     GRINS::AntiochMixtureAveragedTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                                     Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>,
                                                                                                     Antioch::KineticsTheoryThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real>,libMesh::Real>,
                                                                                                     Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> > >(physics_name,input) );
#else
                libmesh_error_msg("ERROR: Antioch requires GSL in order to use kinetics theory based models!");
#endif // ANTIOCH_HAVE_GSL
              }
            else
              this->grins_antioch_model_error_msg(viscosity_model,conductivity_model,diffusivity_model,thermo_model);
          }

        else if( transport_model == std::string("constant") )
          {
            if( viscosity_model != std::string("constant") )
              libmesh_error_msg("Error: For constant transport_model, viscosity model must be constant!");

            if( diffusivity_model != std::string("constant_lewis") )
              libmesh_error_msg("Error: For constant transport_model, diffusivity model must be constant_lewis!");

            if( (thermo_model == std::string("stat_mech")) &&
                (conductivity_model == std::string("constant")) )
              {
                new_physics.reset(new DerivedPhysics<GRINS::AntiochConstantTransportMixture<GRINS::ConstantConductivity>,
                                                     GRINS::AntiochConstantTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>, GRINS::ConstantConductivity> >
                                  (physics_name,input) );
              }
            else if( (thermo_model == std::string("cea")) &&
                     (conductivity_model == std::string("constant")) )
              {
                new_physics.reset(new DerivedPhysics<GRINS::AntiochConstantTransportMixture<GRINS::ConstantConductivity>,
                                                     GRINS::AntiochConstantTransportEvaluator<Antioch::CEAEvaluator<libMesh::Real>, GRINS::ConstantConductivity> >
                                  (physics_name,input) );
              }
            else if( (thermo_model == std::string("stat_mech")) &&
                (conductivity_model == std::string("constant_prandtl")) )
              {
                new_physics.reset(new DerivedPhysics<GRINS::AntiochConstantTransportMixture<GRINS::ConstantPrandtlConductivity>,
                                                     GRINS::AntiochConstantTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>, GRINS::ConstantPrandtlConductivity> >
                                  (physics_name,input) );
              }
            else if( (thermo_model == std::string("cea")) &&
                     (conductivity_model == std::string("constant_prandtl")) )
              {
                new_physics.reset(new DerivedPhysics<GRINS::AntiochConstantTransportMixture<GRINS::ConstantPrandtlConductivity>,
                                                     GRINS::AntiochConstantTransportEvaluator<Antioch::CEAEvaluator<libMesh::Real>, GRINS::ConstantPrandtlConductivity> >
                                  (physics_name,input) );
              }
            else
              this->grins_antioch_model_error_msg(viscosity_model,conductivity_model,diffusivity_model,thermo_model);
          }
        else // transport_model
          {
            std::string error = "Error: Unknown Antioch transport_model "+transport_model+"!\n";
            error += "       Valid values are: constant\n";
            error += "                         mixture_averaged\n";

            libmesh_error_msg(error);
          }
#else
        libmesh_error_msg("Error: Antioch not enabled in this configuration. Reconfigure using --with-antioch option.");

#endif // GRINS_HAVE_ANTIOCH
      }

    else
      {
        std::string error = "Error: Invalid thermo-chemistry library"+thermochem_lib+"!\n";
        error += "       Valid values are: antioch\n";
        error += "                         cantera\n";

        libmesh_error_msg(error);
      }

    libmesh_assert(new_physics);

    return new_physics;
  }

  template<template<typename,typename> class DerivedPhysics>
  void PhysicsFactoryReactingFlows<DerivedPhysics>::grins_antioch_model_error_msg
  ( const std::string& viscosity_model,
    const std::string& conductivity_model,
    const std::string& diffusivity_model,
    const std::string& thermo_model ) const
  {
    std::string error = "Error: Unknown Antioch model combination:\n";
    error += "viscosity_model    = "+viscosity_model+"\n";
    error += "conductivity_model = "+conductivity_model+"\n";
    error += "diffusivity_model = "+diffusivity_model+"\n";
    error += "thermo_model = "+thermo_model+"\n";

    libmesh_error_msg(error);
  }

  ReactingFlowsPhysicsFactoryInitializer::ReactingFlowsPhysicsFactoryInitializer()
  {
    static PhysicsFactoryReactingFlows<ReactingLowMachNavierStokes>
      grins_factory_rlmns
      (PhysicsNaming::reacting_low_mach_navier_stokes(),
       PhysicsNaming::reacting_low_mach_navier_stokes());

    static PhysicsFactoryReactingFlows<ReactingLowMachNavierStokesSPGSMStabilization>
      grins_factory_rlmns_spgsm_stab
      (PhysicsNaming::reacting_low_mach_navier_stokes_spgsm_stab(),
       PhysicsNaming::reacting_low_mach_navier_stokes());
  }

} // end namespace GRINS
