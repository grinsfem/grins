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
#include "grins/materials_parsing.h"

// GRINS
#include "grins/common.h"
#include "grins/grins_physics_names.h"
#include "grins/parameter_user.h"

namespace GRINS
{
  void MaterialsParsing::viscosity_model( const GetPot& input,
                                          const std::string& /*physics*/,
                                          const std::string& material,
                                          std::string& model )
  {
    if( !input.have_variable("Materials/"+material+"/Viscosity/model") )
      {
        libmesh_error_msg("Error: Could not find Materials/"+material+"/Viscosity/model in input file.");
      }

    model = input( "Materials/"+material+"/Viscosity/model", "DIE!" );

    return;
  }

  void MaterialsParsing::thermal_conductivity_model( const GetPot& input,
                                                     const std::string& /*physics*/,
                                                     const std::string& material,
                                                     std::string& model )
  {
    if( !input.have_variable("Materials/"+material+"/ThermalConductivity/model") )
      {
        libmesh_error_msg("Error: Could not find Materials/"+material+"/ThermalConductivity/model in input file.");
      }

    model = input( "Materials/"+material+"/ThermalConductivity/model", "DIE!" );

    return;
  }

  void MaterialsParsing::turb_viscosity_model( const GetPot& input,
                                               const std::string& /*physics*/,
                                               const std::string& material,
                                               std::string& model )
  {
    if( !input.have_variable("Materials/"+material+"/Viscosity/turb_visc_model") )
      {
        libmesh_error_msg("Error: Could not find Materials/"+material+"/Viscosity/turb_visc_model in input file.");
      }

    model = input( "Materials/"+material+"/Viscosity/turb_visc_model", "DIE!" );
  }

  void MaterialsParsing::read_density( const std::string& core_physics_name,
                                       const GetPot& input,
                                       ParameterUser& params,
                                       libMesh::Real& rho )
  {
    std::string material = input("Physics/"+core_physics_name+"/material", "DIE!");

    // Error if both material/Density and rho are specified
    if( input.have_variable("Physics/"+core_physics_name+"/rho") &&
        input.have_variable("Materials/"+material+"/Density/value") )
      {
        libmesh_error_msg("ERROR: Can't specify both rho and Density!");
      }

    // It's deprecated to have nothing and default to 1.0
    if( !input.have_variable("Physics/"+core_physics_name+"/rho") &&
        ( !input.have_variable("Physics/"+core_physics_name+"/material") ||
          !input.have_variable("Materials/"+material+"/Density/value") ) )
      {
        // For some insane reason, we'd originally tied density specifically
        // to incompressible_navier_stokes, so we'll check for that first.
        if( input.have_variable("Physics/"+incompressible_navier_stokes+"/rho") )
          {
            std::string warning = "WARNING: neither Physics/"+core_physics_name+"/rho nor\n";
            warning += "         Physics/"+core_physics_name+"/material options were detected.\n";
            warning += "         But we found Physics/"+incompressible_navier_stokes+"/rho so we using that.\n";
            warning += "        This behavior is DEPRECATED.\n";
            warning += "         Please update and use Physics/"+core_physics_name+"/material.\n";
            grins_warning(warning);

            params.set_parameter
              (rho, input,
               "Physics/"+incompressible_navier_stokes+"/rho", 1.0 /*default*/);
          }
        // Otherwise, the insanity continued and we defaulted to 1.0
        else
          {
            std::string warning = "WARNING: neither Physics/"+core_physics_name+"/rho nor\n";
            warning += "         Physics/"+core_physics_name+"/material  nor\n";
            warning += "         Physics/"+incompressible_navier_stokes+"/rho options were detected.\n";
            warning += "         We are assuming a density value of 1.0. This is DEPRECATED.\n";
            warning += "         Please update and use Physics/"+core_physics_name+"/material.\n";
            grins_warning(warning);

            params.set_parameter
              (rho, input,
               "Physics/"+core_physics_name+"/rho", 1.0 /*default*/);
          }
      }

    // It's deprecated to use rho as the density input
    if( input.have_variable("Physics/"+core_physics_name+"/rho") )
      {
        std::string warning = "WARNING: Using input option Physics/"+core_physics_name+"/rho is DEPRECATED.\n";
        warning += "         Please update and use Physics/"+core_physics_name+"/material.\n";
        grins_warning(warning);

        params.set_parameter
          (rho, input,
           "Physics/"+core_physics_name+"/rho", 1.0 /*default*/);
      }

    // This is the preferred version
    if( input.have_variable("Physics/"+core_physics_name+"/material") &&
        input.have_variable("Materials/"+material+"/Density/value") )
      {
        params.set_parameter
          (rho, input,
           "Materials/"+material+"/Density/value", 0.0 /*default*/);
      }

    // Let's make sure we actually got a valid density value
    if( rho <= 0.0 )
      {
        libmesh_error_msg("ERROR: Detected non-positive value of density!");
      }
  }

} // end namespace GRINS
