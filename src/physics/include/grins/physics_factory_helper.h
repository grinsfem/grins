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

#ifndef GRINS_PHYSICS_FACTORY_HELPER_H
#define GRINS_PHYSICS_FACTORY_HELPER_H

// C++
#include <string>

// libMesh
class GetPot;

namespace GRINS
{

  //! Helper functions for PhysicsFactory
  class PhysicsFactoryHelper
  {
  public:
    PhysicsFactoryHelper(){};
    ~PhysicsFactoryHelper(){};

    //! Determine viscosity model used by turblence classes
    static void parse_turb_viscosity_model( const GetPot& input,
                                            const std::string& physics,
                                            std::string& model );

    //! Determine stress-strain law used by solid mechanics classes
    static void parse_stress_strain_model( const GetPot& input,
                                           const std::string& physics,
                                           std::string& model,
                                           std::string& strain_energy );

    //! Determine thermochemistry model type
    static void parse_thermochemistry_model( const GetPot& input,
                                             const std::string& physics,
                                             std::string& model );

    static void parse_antioch_models( const GetPot& input,
                                      const std::string& physics,
                                      std::string& transport_model,
                                      std::string& thermo_model,
                                      std::string& viscosity_model,
                                      std::string& conductivity_model,
                                      std::string& diffusivity_model );

  private:

    static void deprecated_visc_model_parsing( bool have_viscosity_model,
                                               bool have_material,
                                               const GetPot& input,
                                               const std::string& physics,
                                               std::string& model );

  };

} // end namespace GRINS

#endif // GRINS_PHYSICS_FACTORY_HELPER_H
