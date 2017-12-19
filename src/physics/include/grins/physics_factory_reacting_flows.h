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

#ifndef GRINS_PHYSICS_FACTORY_REACTING_FLOWS_H
#define GRINS_PHYSICS_FACTORY_REACTING_FLOWS_H

// GRINS
#include "grins/physics_factory_with_core.h"
#include "grins/antioch_mixture_averaged_transport_mixture.h"
#include "grins/antioch_mixture_averaged_transport_evaluator.h"
#include "grins/antioch_constant_transport_mixture.h"
#include "grins/antioch_constant_transport_evaluator.h"
#include "grins/antioch_options_naming.h"
#include "grins/antioch_constant_transport_mixture_builder.h"
#include "grins/antioch_mixture_averaged_transport_mixture_builder.h"

namespace GRINS
{
  template<template<typename,typename> class DerivedPhysics>
  class PhysicsFactoryReactingFlows : public PhysicsFactoryWithCore
  {
  public:
    PhysicsFactoryReactingFlows( const std::string& physics_name,
                                 const std::string& core_physics_name )
      : PhysicsFactoryWithCore(physics_name,core_physics_name)
    {}

    ~PhysicsFactoryReactingFlows(){};

  protected:

    virtual std::unique_ptr<Physics> build_physics( const GetPot& input,
                                                    const std::string& physics_name );

    void grins_antioch_model_error_msg( const std::string& viscosity_model,
                                        const std::string& conductivity_model,
                                        const std::string& diffusivity_model ) const;

  private:

    void build_mix_avged_physics( const GetPot & input, const std::string & physics_name,
                                  const std::string & material,
                                  const std::string & thermo_model, const std::string & diffusivity_model,
                                  const std::string & conductivity_model, const std::string & viscosity_model,
                                  std::unique_ptr<Physics> & new_physics );

    void build_const_physics( const GetPot & input, const std::string & physics_name,
                              const std::string & material,
                              const std::string & thermo_model, const std::string & diffusivity_model,
                              const std::string & conductivity_model, const std::string & viscosity_model,
                              std::unique_ptr<Physics> & new_physics );

#ifdef GRINS_HAVE_ANTIOCH
    template<typename KineticsThermo,typename Thermo>
    void build_mix_avged_physics_with_thermo( const GetPot & input, const std::string & physics_name,
                                              const std::string & material,
                                              const std::string & diffusivity_model,
                                              const std::string & conductivity_model,
                                              const std::string & viscosity_model,
                                              std::unique_ptr<Physics> & new_physics )
    {
      if( (diffusivity_model == AntiochOptions::constant_lewis_diffusivity_model()) &&
          (conductivity_model == AntiochOptions::eucken_conductivity_model()) &&
          (viscosity_model == AntiochOptions::sutherland_viscosity_model()) )
        {
          this->build_mix_avged_physics_ptr<KineticsThermo,
                                            Thermo,
                                            Antioch::SutherlandViscosity<libMesh::Real>,
                                            Antioch::EuckenThermalConductivity<Thermo>,
                                            Antioch::ConstantLewisDiffusivity<libMesh::Real> >(input,physics_name,material,new_physics);
        }
      else if( (diffusivity_model == AntiochOptions::constant_lewis_diffusivity_model()) &&
               (conductivity_model == AntiochOptions::eucken_conductivity_model()) &&
               (viscosity_model == AntiochOptions::blottner_viscosity_model()) )
        {
          this->build_mix_avged_physics_ptr<KineticsThermo,
                                            Thermo,
                                            Antioch::BlottnerViscosity<libMesh::Real>,
                                            Antioch::EuckenThermalConductivity<Thermo>,
                                            Antioch::ConstantLewisDiffusivity<libMesh::Real> >(input,physics_name,material,new_physics);
        }
      else if( (diffusivity_model == AntiochOptions::kinetic_theory_diffusivity_model()) &&
               (conductivity_model == AntiochOptions::kinetic_theory_conductivity_model()) &&
               (viscosity_model == AntiochOptions::kinetic_theory_viscosity_model()) )
        {
#ifdef ANTIOCH_HAVE_GSL
          this->build_mix_avged_physics_ptr<KineticsThermo,
                                            Thermo,
                                            Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>,
                                            Antioch::KineticsTheoryThermalConductivity<Thermo,libMesh::Real>,
                                            Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> >(input,physics_name,material,new_physics);
#else
          libmesh_error_msg("ERROR: Antioch requires GSL in order to use kinetics theory based models!");
#endif // ANTIOCH_HAVE_GSL
        }
      else
        this->grins_antioch_model_error_msg(viscosity_model,conductivity_model,diffusivity_model);
    }

    template<typename KineticsThermo,typename Thermo>
    void build_const_physics_with_thermo( const GetPot & input, const std::string & physics_name,
                                          const std::string & material, const std::string & conductivity_model,
                                          std::unique_ptr<Physics> & new_physics )
    {
      if( conductivity_model == AntiochOptions::constant_conductivity_model() )
        {
          this->build_const_physics_ptr<KineticsThermo,Thermo,ConstantConductivity>
            (input,physics_name,material,new_physics);
        }
      else if( conductivity_model == AntiochOptions::constant_prandtl_conductivity_model() )
        {
          this->build_const_physics_ptr<KineticsThermo,Thermo,ConstantPrandtlConductivity>
            (input,physics_name,material,new_physics);
        }
      else
        {
          std::string error = "ERROR: Invalid conductivity model for constant transport!\n";
          error += "       Valid choices are "+AntiochOptions::constant_conductivity_model()+"\n";
          error += "                         "+AntiochOptions::constant_prandtl_conductivity_model()+"\n";

          libmesh_error_msg(error);
        }
    }

    template<typename KineticsThermo,typename Thermo,typename Viscosity,typename Conductivity,typename Diffusivity>
    void build_mix_avged_physics_ptr( const GetPot& input, const std::string& physics_name,
                                      const std::string & material, std::unique_ptr<Physics> & new_physics )
    {
      AntiochMixtureAveragedTransportMixtureBuilder mix_builder;

      std::unique_ptr<AntiochMixtureAveragedTransportMixture<KineticsThermo,Thermo,Viscosity,Conductivity,Diffusivity> >
        gas_mixture = mix_builder.build_mixture<KineticsThermo,Thermo,Viscosity,Conductivity,Diffusivity>(input,material);

      new_physics.reset(new DerivedPhysics<AntiochMixtureAveragedTransportMixture<KineticsThermo,Thermo,Viscosity,Conductivity,Diffusivity>,
                        AntiochMixtureAveragedTransportEvaluator<KineticsThermo,Thermo,Viscosity,Conductivity,Diffusivity> >
                        (physics_name,input,gas_mixture) );
    }

    template<typename KineticsThermo,typename Thermo,typename Conductivity>
    void build_const_physics_ptr( const GetPot& input, const std::string& physics_name,
                                  const std::string & material, std::unique_ptr<Physics> & new_physics )
    {
      AntiochConstantTransportMixtureBuilder mix_builder;

      std::unique_ptr<GRINS::AntiochConstantTransportMixture<KineticsThermo,Conductivity> >
        gas_mixture = mix_builder.build_mixture<KineticsThermo,Conductivity>(input,material);

      new_physics.reset(new DerivedPhysics<AntiochConstantTransportMixture<KineticsThermo,Conductivity>,
                        AntiochConstantTransportEvaluator<KineticsThermo,Thermo,Conductivity> >
                        (physics_name,input,gas_mixture) );
    }
#endif // GRINS_HAVE_ANTIOCH

  };

  class ReactingFlowsPhysicsFactoryInitializer
  {
  public:
    ReactingFlowsPhysicsFactoryInitializer();
    ~ReactingFlowsPhysicsFactoryInitializer(){}
  };

} // end namespace GRINS

#endif // GRINS_PHYSICS_FACTORY_REACTING_FLOWS_H
