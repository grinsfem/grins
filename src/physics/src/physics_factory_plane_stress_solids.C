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

// This class
#include "grins/physics_factory_plane_stress_solids.h"

// GRINS
#include "grins/physics_factory_helper.h"
#include "grins/hookes_law.h"
#include "grins/hookes_law_1d.h"
#include "grins/incompressible_plane_stress_hyperelasticity.h"
#include "grins/mooney_rivlin.h"

// Physics whose factories we're instantiating
#include "grins/elastic_membrane.h"
#include "grins/elastic_membrane_rayleigh_damping.h"

namespace GRINS
{
  template<template<typename> class DerivedPhysics>
  libMesh::UniquePtr<Physics> PhysicsFactoryPlaneStressSolids<DerivedPhysics>::build_physics( const GetPot& input,
                                                                                              const std::string& physics_name )
  {
    std::string core_physics = this->find_core_physics_name(physics_name);

    std::string model = "none";
    std::string strain_energy = "none";

    PhysicsFactoryHelper::parse_stress_strain_model( input,
                                                     core_physics,
                                                     model,
                                                     strain_energy );

    libMesh::UniquePtr<Physics> new_physics;

    if( model == std::string("hookes_law") )
      new_physics.reset( new DerivedPhysics<HookesLaw>
                         (physics_name,input, false /*is_compressible*/) );

    else if( model == std::string("incompressible_hyperelasticity") )
      {
        if( strain_energy == std::string("mooney_rivlin") )
          {
            new_physics.reset( new DerivedPhysics<IncompressiblePlaneStressHyperelasticity<MooneyRivlin> >(physics_name,input,false /*is_compressible*/) );
          }
        else
          {
            std::string error = "ERROR: Invalid strain_energy "+strain_energy+"!\n";
            error += "       Valid values are: mooney_rivlin\n";
            libmesh_error_msg(error);
          }

      }
    else
      {
        std::string error = "Error: Invalid stress-strain model: "+model+"!\n";
        error += "       Valid values are: hookes_law\n";
        error += "                         incompressible_hyperelasticity\n";
        libmesh_error_msg(error);
      }

    libmesh_assert(new_physics);

    return new_physics;
  }

  // Instantiate all the "variable density flow" Physics factories.
  PhysicsFactoryPlaneStressSolids<ElasticMembrane> grins_factory_elastic_membrane
  (PhysicsNaming::elastic_membrane(),PhysicsNaming::elastic_membrane());

  PhysicsFactoryPlaneStressSolids<ElasticMembraneRayleighDamping> grins_factory_elastic_membrane_rayleigh_damping
  (PhysicsNaming::elastic_membrane_rayleigh_damping(),PhysicsNaming::elastic_membrane());

} // end namespace GRINS
