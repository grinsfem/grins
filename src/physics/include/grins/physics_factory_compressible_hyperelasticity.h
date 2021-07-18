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

#ifndef GRINS_PHYSICS_FACTORY_COMPRESSIBLE_HYPERELASTICITY_H
#define GRINS_PHYSICS_FACTORY_COMPRESSIBLE_HYPERELASTICITY_H

// GRINS
#include "grins/materials_parsing.h"
#include "grins/physics_factory_with_core.h"
#include "grins/compressible_mooney_rivlin.h"

namespace GRINS
{
  template<unsigned int Dim,template<unsigned int,typename> class DerivedPhysics>
  class PhysicsFactoryCompressibleHyperelasticity : public PhysicsFactoryWithCore
  {
  public:
    PhysicsFactoryCompressibleHyperelasticity( const std::string& physics_name,
                                               const std::string& core_physics_name )
      : PhysicsFactoryWithCore(physics_name,core_physics_name)
    {}

    ~PhysicsFactoryCompressibleHyperelasticity(){};

  protected:

    virtual std::unique_ptr<Physics> build_physics( const GetPot& input,
                                                    const std::string& physics_name ) override;

  };

  template<unsigned int Dim,template<unsigned int,typename> class DerivedPhysics>
  inline
  std::unique_ptr<Physics>
  PhysicsFactoryCompressibleHyperelasticity<Dim,DerivedPhysics>::build_physics
  ( const GetPot& input, const std::string& physics_name )
  {
    std::string core_physics = this->find_core_physics_name(physics_name);

    std::string strain_energy("none");
    MaterialsParsing::strain_energy(input,physics_name,strain_energy);

    std::unique_ptr<Physics> new_physics;

    if( strain_energy == std::string("compressible_mooney_rivlin") )
      new_physics.reset( new DerivedPhysics<Dim,CompressibleMooneyRivlin>(physics_name,input) );
    else
      {
        std::string error = "ERROR: Invalid strain_energy "+strain_energy+"!\n";
        error += "       Valid values are: compressible_mooney_rivlin\n";
        libmesh_error_msg(error);
      }

    libmesh_assert(new_physics);

    return new_physics;
  }

} // end namespace GRINS

#endif //GRINS_PHYSICS_FACTORY_COMPRESSIBLE_HYPERELASTICITY_H
