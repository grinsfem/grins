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

} // end namespace GRINS