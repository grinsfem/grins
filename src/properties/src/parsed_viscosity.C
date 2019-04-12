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


// This class
#include "grins/parsed_viscosity.h"

//GRINS
#include "grins/physics_naming.h"
#include "grins/common.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{

  ParsedViscosity::ParsedViscosity( const GetPot& input )
    : ParsedPropertyBase(),
      ParameterUser("ParsedViscosity"),
      ViscosityBase()
  {

    // Warning about this constructor being deprecated
    {
      std::string warning = "WARNING: Use of this constructor is DEPRECATED.\n";
      warning += "         Please update to use constructor with input material name.\n";
      grins_warning(warning);
    }

    this->set_parameter(this->_func, input,
                        "Materials/Viscosity/mu",
                        "DIE!");

    std::string viscosity_function = input("Materials/Viscosity/mu",std::string("0"));

    if( !this->check_func_nonzero(viscosity_function) )
      {
        libmesh_error_msg("ERROR: Detected '0' function for ParsedConductivity!");
      }
  }

  ParsedViscosity::ParsedViscosity( const GetPot& input, const std::string& material )
    : ParameterUser("ParsedViscosity"),
      ViscosityBase()
  {
    this->check_input_consistency(input,material);

    std::string viscosity_function = "0";

    // If we have the new version, we parse that
    if( input.have_variable("Materials/"+material+"/Viscosity/value") )
      {
        this->set_parameter(this->_func, input,
                            "Materials/"+material+"/Viscosity/value",
                            "DIE!");

        viscosity_function = input("Materials/"+material+"/Viscosity/value",std::string("0"));
      }
    // If we have the old DEPRECATED version, use that
    else if( input.have_variable("Materials/Viscosity/mu") )
      {
        this->old_mu_warning();

        this->set_parameter(this->_func, input,
                            "Materials/Viscosity/mu",
                            "DIE!");

        viscosity_function = input("Materials/Viscosity/mu",std::string("0"));
      }
    // If we don't have either, that's an error
    else
      {
        libmesh_error_msg("Error: Could not find either Materials/"+material+"/Viscosity/value or Materials/Viscosity/mu");
      }

    if( !this->check_func_nonzero(viscosity_function) )
      {
        libmesh_error_msg("ERROR: Detected '0' function for ParsedConductivity!");
      }
  }

  ParsedViscosity::~ParsedViscosity()
  {
    return;
  }

} // namespace GRINS
