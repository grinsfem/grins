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

#ifndef GRINS_PHYSICS_FACTORY_VARIABLE_DENSITY_FLOW_H
#define GRINS_PHYSICS_FACTORY_VARIABLE_DENSITY_FLOW_H

// GRINS
#include "grins/materials_parsing.h"
#include "grins/physics_factory_with_core.h"
#include "grins/constant_viscosity.h"
#include "grins/constant_conductivity.h"
#include "grins/constant_specific_heat.h"

namespace GRINS
{
  template<template<typename,typename,typename> class DerivedPhysics>
  class PhysicsFactoryVariableDensityFlow : public PhysicsFactoryWithCore
  {
  public:
    PhysicsFactoryVariableDensityFlow( const std::string& physics_name,
                                       const std::string& core_physics_name )
      : PhysicsFactoryWithCore(physics_name,core_physics_name)
    {}

    ~PhysicsFactoryVariableDensityFlow(){};

  protected:

    virtual std::unique_ptr<Physics> build_physics( const GetPot& input,
                                                    const std::string& physics_name );

    void prop_error_msg( const std::string& physics,
                         const std::string& conductivity,
                         const std::string& viscosity,
                         const std::string& specific_heat ) const;

  };

  template<template<typename,typename,typename> class DerivedPhysics>
  inline
  std::unique_ptr<Physics> PhysicsFactoryVariableDensityFlow<DerivedPhysics>::build_physics( const GetPot& input,
                                                                                             const std::string& physics_name )
  {
    std::string core_physics = this->find_core_physics_name(physics_name);

    std::string conductivity;
    MaterialsParsing::thermal_conductivity_model(input,core_physics,conductivity);

    std::string viscosity;
    MaterialsParsing::viscosity_model(input,core_physics,viscosity);

    std::string specific_heat;
    MaterialsParsing::specific_heat_model(input,core_physics,specific_heat);

    std::unique_ptr<Physics> new_physics;

    if(  conductivity == "constant" && viscosity == "constant" && specific_heat == "constant" )
      new_physics.reset( new DerivedPhysics<ConstantViscosity,ConstantSpecificHeat,ConstantConductivity>
                         (physics_name,input) );

    else
      this->prop_error_msg(physics_name, conductivity, viscosity, specific_heat);

    libmesh_assert(new_physics);

    return new_physics;
  }

  template<template<typename,typename,typename> class DerivedPhysics>
  inline
  void PhysicsFactoryVariableDensityFlow<DerivedPhysics>::prop_error_msg( const std::string& physics,
                                                                          const std::string& conductivity,
                                                                          const std::string& viscosity,
                                                                          const std::string& specific_heat ) const
  {
    std::string error = "================================================================\n";
    error += "Invalid combination of models for "+physics+"\n";
    error += "Viscosity model     = "+viscosity+"\n";
    error += "Conductivity model  = "+conductivity+"\n";
    error += "Specific heat model = "+specific_heat+"\n";
    error += "================================================================\n";

    libmesh_error_msg(error);
  }

} // end namespace GRINS

#endif // GRINS_PHYSICS_FACTORY_VARIABLE_DENSITY_FLOW_H
