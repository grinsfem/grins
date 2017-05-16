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

#ifndef GRINS_PHYSICS_FACTORY_REACTING_FLOWS_H
#define GRINS_PHYSICS_FACTORY_REACTING_FLOWS_H

// GRINS
#include "grins/physics_factory_with_core.h"
#include "grins/antioch_mixture_averaged_transport_mixture.h"
#include "grins/antioch_mixture_averaged_transport_evaluator.h"

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

    virtual libMesh::UniquePtr<Physics> build_physics( const GetPot& input,
                                                       const std::string& physics_name );

    void grins_antioch_model_error_msg( const std::string& viscosity_model,
                                        const std::string& conductivity_model,
                                        const std::string& diffusivity_model ) const;

  private:

    void build_mix_avged_physics( const GetPot & input, const std::string & physics_name,
                                  const std::string & thermo_model, const std::string & diffusivity_model,
                                  const std::string & conductivity_model, const std::string & viscosity_model,
                                  libMesh::UniquePtr<Physics> & new_physics );

    template<typename KineticsThermo,typename Thermo>
    void build_mix_avged_physics_with_thermo( const GetPot & input, const std::string & physics_name,
                                              const std::string & diffusivity_model,
                                              const std::string & conductivity_model,
                                              const std::string & viscosity_model,
                                              libMesh::UniquePtr<Physics> & new_physics )
    {
      if( (diffusivity_model == std::string("constant_lewis")) &&
          (conductivity_model == std::string("eucken")) &&
          (viscosity_model == std::string("sutherland")) )
      {
        this->build_mix_avged_physics_ptr<KineticsThermo,
                                          Thermo,
                                          Antioch::SutherlandViscosity<libMesh::Real>,
                                          Antioch::EuckenThermalConductivity<Thermo>,
                                          Antioch::ConstantLewisDiffusivity<libMesh::Real> >(input,physics_name,new_physics);
      }
    else if( (diffusivity_model == std::string("constant_lewis")) &&
             (conductivity_model == std::string("eucken")) &&
             (viscosity_model == std::string("blottner")) )
      {
        this->build_mix_avged_physics_ptr<KineticsThermo,
                                          Thermo,
                                          Antioch::BlottnerViscosity<libMesh::Real>,
                                          Antioch::EuckenThermalConductivity<Thermo>,
                                          Antioch::ConstantLewisDiffusivity<libMesh::Real> >(input,physics_name,new_physics);
      }
    else if( (diffusivity_model == std::string("kinetics_theory")) &&
             (conductivity_model == std::string("kinetics_theory")) &&
             (viscosity_model == std::string("kinetics_theory")) )
      {
#ifdef ANTIOCH_HAVE_GSL
        this->build_mix_avged_physics_ptr<KineticsThermo,
                                          Thermo,
                                          Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>,
                                          Antioch::KineticsTheoryThermalConductivity<Thermo,libMesh::Real>,
                                          Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> >(input,physics_name,new_physics);
#else
        libmesh_error_msg("ERROR: Antioch requires GSL in order to use kinetics theory based models!");
#endif // ANTIOCH_HAVE_GSL
      }
    else
      this->grins_antioch_model_error_msg(viscosity_model,conductivity_model,diffusivity_model);
    }

    template<typename KineticsThermo,typename Thermo,typename Viscosity,typename Conductivity,typename Diffusivity>
    void build_mix_avged_physics_ptr( const GetPot& input, const std::string& physics_name,
                                      libMesh::UniquePtr<Physics> & new_physics )
    {
      new_physics.reset(new DerivedPhysics<AntiochMixtureAveragedTransportMixture<KineticsThermo,Thermo,Viscosity,Conductivity,Diffusivity>,
                                           AntiochMixtureAveragedTransportEvaluator<KineticsThermo,Thermo,Viscosity,Conductivity,Diffusivity> >
                        (physics_name,input) );
    }

  };

  class ReactingFlowsPhysicsFactoryInitializer
  {
  public:
    ReactingFlowsPhysicsFactoryInitializer();
    ~ReactingFlowsPhysicsFactoryInitializer(){}
  };

} // end namespace GRINS

#endif // GRINS_PHYSICS_FACTORY_REACTING_FLOWS_H
