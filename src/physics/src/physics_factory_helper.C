//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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

// These functions
#include "grins/physics_factory_helper.h"

// GRINS
#include "grins/common.h"
#include "grins/grins_physics_names.h"
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  void PhysicsFactoryHelper::parse_viscosity_model( const GetPot& input,
                                                    const std::string& physics,
                                                    std::string& model )
  {
    // Newer, preferred version
    bool have_material = MaterialsParsing::have_material(input,physics);

    // Old deprecated version
    bool have_viscosity_model = input.have_variable("Physics/"+incompressible_navier_stokes+"/viscosity_model");

    if( have_material && have_viscosity_model )
      {
        libmesh_error_msg("Error: Cannot specify both viscosity_model and material.");
      }

    /* If the user hasn't specified the material or the viscosity_model,
       we're assuming they're using the old version. */
    if( !have_viscosity_model && !have_material )
      {
        std::string warning = "Warning: Neither viscosity_model nor material were specified.\n";
        warning += "      We are assuming a constant viscosity model.\n";
        warning += "      This case is DEPRECATED.\n";
        warning += "      Please update and specify Physics/"+physics+"/material.\n";
        grins_warning(warning);

        model = input( "Physics/"+incompressible_navier_stokes+"/viscosity_model", "constant" );
      }

    if( have_viscosity_model )
      {
        std::string warning = "Warning: Option Physics/"+incompressible_navier_stokes+"/viscosity_model is DEPRECATED.\n";
        warning += "         Please update to use Physics/"+physics+"/material.\n";
        grins_warning(warning);

        model = input( "Physics/"+incompressible_navier_stokes+"/viscosity_model", "constant" );
      }

    if( have_material )
      {
        std::string material;
        MaterialsParsing::material_name( input, physics, material );
        MaterialsParsing::viscosity_model( input, physics, material, model );
      }

    return;
  }

  void PhysicsFactoryHelper::parse_conductivity_model( const GetPot& input,
                                                       const std::string& physics,
                                                       std::string& model )
  {
    // Newer, preferred version
    bool have_material = MaterialsParsing::have_material(input,physics);

    // Old deprecated version
    bool have_conductivity_model = input.have_variable("Physics/"+heat_transfer+"/conductivity_model");

    if( have_material && have_conductivity_model )
      {
        libmesh_error_msg("Error: Cannot specify both conductivity_model and material.");
      }

    /* If the user hasn't specified the material or the conductivity_model,
       we're assuming they're using the old version.
       This is deprecated.*/
    if( !have_conductivity_model && !have_material )
      {
        std::string warning = "Warning: Neither conductivity_model nor material were specified.\n";
        warning += "      We are assuming a constant conductivity model.\n";
        warning += "      This case is DEPRECATED.\n";
        warning += "      Please update and specify Physics/"+physics+"/material.\n";
        grins_warning(warning);

        model = input( "Physics/"+heat_transfer+"/conductivity_model", "constant" );
      }

    // Deprecated
    if( have_conductivity_model )
      {
        std::string warning = "Warning: Option Physics/"+heat_transfer+"/conductivity_model is DEPRECATED.\n";
        warning += "         Please update to use Physics/"+physics+"/material.\n";
        grins_warning(warning);

        model = input( "Physics/"+heat_transfer+"/conductivity_model", "constant" );
      }

    // Preferred
    if( have_material )
      {
        std::string material;
        MaterialsParsing::material_name( input, physics, material );
        MaterialsParsing::thermal_conductivity_model( input, physics, material, model );
      }

    return;
  }

  void PhysicsFactoryHelper::parse_specific_heat_model( const GetPot& input,
                                                        const std::string& physics,
                                                        std::string& model )
  {
    model = input( "Physics/"+low_mach_navier_stokes+"/specific_heat_model", "constant" );
    return;
  }

} // end namespace GRINS
