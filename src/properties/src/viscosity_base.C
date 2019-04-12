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
#include "grins/viscosity_base.h"

// GRINS
#include "grins/common.h"
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  void ViscosityBase::check_input_consistency( const GetPot& input, const std::string& material ) const
  {
    // We can't have both the materials version and the old versions
    MaterialsParsing::duplicate_input_test(input,
                                           "Materials/"+material+"/Viscosity/value",
                                           "Materials/Viscosity/mu" );

    // If the material section exists, but not the variable, this is an error
    if( input.have_section("Materials/"+material+"/Viscosity") &&
        !input.have_variable("Materials/"+material+"/Viscosity/value") )
      {
        libmesh_error_msg("Error: Found section Materials/"+material+"/Viscosity, but not variable value.");
      }

    return;
  }

  void ViscosityBase::old_mu_warning() const
  {
    std::string warning = "WARNING: Use of Materials/Viscosity/mu is DEPRECATED.\n";
    warning += "         Please update to use  Materials/MATERIAL_NAME/Viscosity/value,\n";
    warning += "         where MATERIAL_NAME is given by Physics/PHYSICS_CLASS/material.\n";
    grins_warning(warning);

    return;
  }

} // end namespace GRINS
