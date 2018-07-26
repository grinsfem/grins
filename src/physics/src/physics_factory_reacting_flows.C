// This class
#include "grins/physics_factory_reacting_flows.h"

// GRINS
#include "grins/common.h"
#include "grins/materials_parsing.h"
#include "grins/cantera_mixture.h"
#include "grins/cantera_evaluator.h"

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
  std::unique_ptr<Physics> PhysicsFactoryReactingFlows<DerivedPhysics>::build_physics( const GetPot& input,
                                                                                       const std::string& physics_name )
  {
    std::string core_physics = this->find_core_physics_name(physics_name);

    std::string material = MaterialsParsing::material_name(input,core_physics);

    std::string thermochem_lib;
    MaterialsParsing::thermochemistry_lib( input, core_physics, thermochem_lib );

    std::unique_ptr<Physics> new_physics;

    if( thermochem_lib == "cantera" )
      {
#ifdef GRINS_HAVE_CANTERA
        std::unique_ptr<CanteraMixture> gas_mix( new CanteraMixture(input,material) );

        new_physics.reset(new DerivedPhysics<CanteraMixture,CanteraEvaluator>(physics_name,input,gas_mix));
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

        MaterialsParsing::antioch_models( input,
                                          core_physics,
                                          transport_model,
                                          thermo_model,
                                          viscosity_model,
                                          conductivity_model,
                                          diffusivity_model );

        if( transport_model == AntiochOptions::mix_avged_transport_model() )
          {
            this->build_mix_avged_physics( input, physics_name, material, thermo_model, diffusivity_model,
                                           conductivity_model, viscosity_model, new_physics );
          }

        else if( transport_model == AntiochOptions::constant_transport_model() )
          {
            this->build_const_physics( input, physics_name, material, thermo_model, diffusivity_model,
                                       conductivity_model, viscosity_model, new_physics );
          }
        else // transport_model
          {
            std::string error = "Error: Unknown Antioch transport_model "+transport_model+"!\n";
            error += "       Valid values are: "+AntiochOptions::constant_transport_model()+"\n";
            error += "                         "+AntiochOptions::mix_avged_transport_model()+"\n";

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



#ifdef GRINS_HAVE_ANTIOCH
  template<template<typename,typename> class DerivedPhysics>
  void PhysicsFactoryReactingFlows<DerivedPhysics>::
  build_mix_avged_physics( const GetPot & input, const std::string & physics_name, const std::string & material,
                           const std::string & thermo_model, const std::string & diffusivity_model,
                           const std::string & conductivity_model, const std::string & viscosity_model,
                           std::unique_ptr<Physics> & new_physics )
  {
    AntiochMixtureBuilderBase builder;

    if( (thermo_model == AntiochOptions::stat_mech_thermo_model()) )
      {
        ThermoEnum thermo_type = builder.get_thermo_type(input,material);
        switch(thermo_type)
          {
          case(NASA7):
            {
              this->build_mix_avged_physics_with_thermo<Antioch::NASA7CurveFit<libMesh::Real>,
                                                        Antioch::StatMechThermodynamics<libMesh::Real> >
                (input,physics_name,material,diffusivity_model,
                 conductivity_model,viscosity_model,
                 new_physics);

              break;
            }
          case(NASA9):
            {
              this->build_mix_avged_physics_with_thermo<Antioch::NASA9CurveFit<libMesh::Real>,
                                                        Antioch::StatMechThermodynamics<libMesh::Real> >
                (input,physics_name,material,diffusivity_model,
                 conductivity_model,viscosity_model,
                 new_physics);

              break;
            }
          case(CEA):
            {
              this->build_mix_avged_physics_with_thermo<Antioch::CEACurveFit<libMesh::Real>,
                                                        Antioch::StatMechThermodynamics<libMesh::Real> >
                (input,physics_name,material,diffusivity_model,
                 conductivity_model,viscosity_model,
                 new_physics);

              break;
            }
          case(INVALID):
          default:
            libmesh_error_msg("ERROR: Invalid thermo type for thermo_model!");
          }
      }
    else if( thermo_model == AntiochOptions::ideal_gas_thermo_model() )
      {
        ThermoEnum thermo_type = builder.get_thermo_type(input,material);
        switch(thermo_type)
          {
          case(NASA7):
            {
              this->build_mix_avged_physics_with_thermo<Antioch::NASA7CurveFit<libMesh::Real>,
                                                        Antioch::IdealGasThermo<Antioch::NASA7CurveFit<libMesh::Real>,libMesh::Real> >
                (input,physics_name,material,diffusivity_model,
                 conductivity_model,viscosity_model,
                 new_physics);

              break;
            }
          case(NASA9):
            {
              this->build_mix_avged_physics_with_thermo<Antioch::NASA9CurveFit<libMesh::Real>,
                                                        Antioch::IdealGasThermo<Antioch::NASA9CurveFit<libMesh::Real>,libMesh::Real> >
                (input,physics_name,material,diffusivity_model,
                 conductivity_model,viscosity_model,
                 new_physics);

              break;
            }
          case(CEA):
            {
              this->build_mix_avged_physics_with_thermo<Antioch::CEACurveFit<libMesh::Real>,
                                                        Antioch::IdealGasThermo<Antioch::CEACurveFit<libMesh::Real>,libMesh::Real> >
                (input,physics_name,material,diffusivity_model,
                 conductivity_model,viscosity_model,
                 new_physics);

              break;
            }
          case(INVALID):
          default:
            libmesh_error_msg("ERROR: Invalid thermo type for thermo_model!");
          }
      }
    else if( (thermo_model == AntiochOptions::cea_nasa_model()) )
      {
        {
          std::string msg = "WARNING: Specifying thermo_model = "+AntiochOptions::cea_nasa_model();
          msg += "is DEPREACATED!\n";
          msg += "         You should specify either "+AntiochOptions::stat_mech_thermo_model();
          msg += "or "+AntiochOptions::ideal_gas_thermo_model()+"\n";
          grins_warning(msg);
        }

        this->build_mix_avged_physics_with_thermo<Antioch::CEACurveFit<libMesh::Real>,
                                                  Antioch::IdealGasThermo<Antioch::CEACurveFit<libMesh::Real>,libMesh::Real> >
          (input,physics_name,material,diffusivity_model,
           conductivity_model,viscosity_model,
           new_physics);
      }
    else
      {
        std::string error = "Error: Unknown Antioch thermo model "+thermo_model+"\n";
        error += "       Valid values are: "+AntiochOptions::stat_mech_thermo_model()+"\n";
        error += "       Valid values are: "+AntiochOptions::cea_nasa_model()+"\n";
      }
  }

  template<template<typename,typename> class DerivedPhysics>
  void PhysicsFactoryReactingFlows<DerivedPhysics>::
  build_const_physics( const GetPot & input, const std::string & physics_name, const std::string & material,
                       const std::string & thermo_model, const std::string & diffusivity_model,
                       const std::string & conductivity_model, const std::string & viscosity_model,
                       std::unique_ptr<Physics> & new_physics )
  {
    // First check the things we must have for constant transport models
    if( viscosity_model != AntiochOptions::constant_viscosity_model() )
      libmesh_error_msg("Error: For constant transport_model, viscosity model must be constant!");

    if( diffusivity_model != AntiochOptions::constant_lewis_diffusivity_model() )
      libmesh_error_msg("Error: For constant transport_model, diffusivity model must be constant_lewis!");


    AntiochMixtureBuilderBase builder;
    ThermoEnum thermo_type = builder.get_thermo_type(input,material);

    if( thermo_model == AntiochOptions::stat_mech_thermo_model() )
      {
        switch(thermo_type)
          {
          case(NASA7):
            {
              this->build_const_physics_with_thermo<Antioch::NASA7CurveFit<libMesh::Real>,
                                                    Antioch::StatMechThermodynamics<libMesh::Real> >
                (input,physics_name,material,conductivity_model,new_physics);

              break;
            }
          case(NASA9):
            {
              this->build_const_physics_with_thermo<Antioch::NASA9CurveFit<libMesh::Real>,
                                                    Antioch::StatMechThermodynamics<libMesh::Real> >
                (input,physics_name,material,conductivity_model,new_physics);

              break;
            }
          case(CEA):
            {
              this->build_const_physics_with_thermo<Antioch::CEACurveFit<libMesh::Real>,
                                                    Antioch::StatMechThermodynamics<libMesh::Real> >
                (input,physics_name,material,conductivity_model,new_physics);

              break;
            }
          case(INVALID):
          default:
            libmesh_error_msg("ERROR: Invalid thermo type for thermo_model!");
          }
      }
    else if( thermo_model == AntiochOptions::ideal_gas_thermo_model() )
      {
        switch(thermo_type)
          {
          case(NASA7):
            {
              this->build_const_physics_with_thermo<Antioch::NASA7CurveFit<libMesh::Real>,
                                                    Antioch::IdealGasThermo<Antioch::NASA7CurveFit<libMesh::Real> > >
                (input,physics_name,material,conductivity_model,new_physics);

              break;
            }
          case(NASA9):
            {
              this->build_const_physics_with_thermo<Antioch::NASA9CurveFit<libMesh::Real>,
                                                    Antioch::IdealGasThermo<Antioch::NASA9CurveFit<libMesh::Real> > >
                (input,physics_name,material,conductivity_model,new_physics);

              break;
            }
          case(CEA):
            {
              this->build_const_physics_with_thermo<Antioch::CEACurveFit<libMesh::Real>,
                                                    Antioch::IdealGasThermo<Antioch::CEACurveFit<libMesh::Real> > >
                (input,physics_name,material,conductivity_model,new_physics);

              break;
            }
          case(INVALID):
          default:
            libmesh_error_msg("ERROR: Invalid thermo type for thermo_model!");
          }
      }
    else if( thermo_model == AntiochOptions::cea_nasa_model() )
      {
        {
          std::string msg = "WARNING: Specifying thermo_model = "+AntiochOptions::cea_nasa_model();
          msg += "is DEPREACATED!\n";
          msg += "         You should specify either "+AntiochOptions::stat_mech_thermo_model();
          msg += "or "+AntiochOptions::ideal_gas_thermo_model()+"\n";
          grins_warning(msg);
        }

        this->build_const_physics_with_thermo<Antioch::CEACurveFit<libMesh::Real>,
                                              Antioch::IdealGasThermo<Antioch::CEACurveFit<libMesh::Real>,libMesh::Real> >
          (input,physics_name,material,conductivity_model,new_physics);
      }
    else
      this->grins_antioch_model_error_msg(viscosity_model,conductivity_model,diffusivity_model);
  }
#endif // GRINS_HAVE_ANTIOCH

  template<template<typename,typename> class DerivedPhysics>
  void PhysicsFactoryReactingFlows<DerivedPhysics>::grins_antioch_model_error_msg
  ( const std::string& viscosity_model,
    const std::string& conductivity_model,
    const std::string& diffusivity_model ) const
  {
    std::string error = "Error: Unknown Antioch model combination:\n";
    error += "viscosity_model    = "+viscosity_model+"\n";
    error += "conductivity_model = "+conductivity_model+"\n";
    error += "diffusivity_model = "+diffusivity_model+"\n";

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
