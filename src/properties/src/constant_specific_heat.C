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
#include "grins/constant_specific_heat.h"

// GRINS
#include "grins/common.h"
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{

  ConstantSpecificHeat::ConstantSpecificHeat( const GetPot& input )
    : ParameterUser("ConstantSpecificHeat"),
      _cp(0.0)
  {
    // Warning about this constructor being deprecated
    {
      std::string warning = "WARNING: Use of this constructor is DEPRECATED.\n";
      warning += "         Please update to use constructor with input material name.\n";
      grins_warning(warning);
    }

    if( !input.have_variable("Materials/SpecificHeat/cp") )
      {
        std::cerr << "Error: Must specify cp value for constant specific heat model!" << std::endl;
        libmesh_error();
      }

    this->set_parameter
      (_cp, input, "Materials/SpecificHeat/cp", _cp);
  }

  ConstantSpecificHeat::ConstantSpecificHeat( const GetPot& input,
                                              const std::string& material )
    : ParameterUser("ConstantSpecificHeat"),
      _cp(0.0)
  {
    // It's an error to have both the old and new version
    MaterialsParsing::duplicate_input_test(input,
                                           "Materials/"+material+"/SpecificHeat/value",
                                           "Materials/SpecificHeat/cp");

    // If we have the "new" version, then parse it
    if( input.have_variable("Materials/"+material+"/SpecificHeat/value") )
      {
        this->set_parameter
          (_cp, input, "Materials/"+material+"/SpecificHeat/value", _cp);
      }
    // If instead we have the old version, use that.
    else if( input.have_variable("Materials/SpecificHeat/cp") )
      {
        MaterialsParsing::dep_input_warning( "Materials/SpecificHeat/cp",
                                             "SpecificHeat/value" );

        this->set_parameter
          (_cp, input, "Materials/SpecificHeat/cp", _cp);
      }
    else
      {
        libmesh_error_msg("ERROR: Could not find valid input for ConstantSpecificHeat! Please set Materials/"+material+"/SpecificHeat/value");
      }
  }

  ConstantSpecificHeat::~ConstantSpecificHeat()
  {
    return;
  }

} // namespace GRINS
