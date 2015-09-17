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


// This class
#include "grins/constant_conductivity.h"

//GRINS
#include "grins/common.h"
#include "grins/grins_physics_names.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{

  ConstantConductivity::ConstantConductivity( const GetPot& input )
    : ParameterUser("ConstantConductivity"),
      _k(0.0)
  {
    // Warning about this constructor being deprecated
    {
      std::string warning = "WARNING: Use of this constructor is DEPRECATED.\n";
      warning += "         Please update to use constructor with input material name.\n";
      grins_warning(warning);
    }

    if( !input.have_variable("Materials/Conductivity/k") )
      {
        libmesh_warning("No Materials/Conductivity/k specified!\n");

	// Try and get the conductivity from other specifications
        this->set_parameter
	  (_k, input, "Physics/"+incompressible_navier_stokes+"/k", _k);
      }
    else
      this->set_parameter
        (_k, input, "Materials/Conductivity/k", _k);
    return;
  }

  ConstantConductivity::ConstantConductivity( const GetPot& input, const std::string& material )
  : ParameterUser("ConstantConductivity"),
    _k(0.0)
  {
    if( input.have_variable("Materials/"+material+"/ThermalConductivity/value") &&
        input.have_variable("Materials/Conductivity/k") )
      {
        libmesh_error_msg("Error: Cannot specify both Materials/"+material+"/ThermalConductivity/value and Materials/Conductivity/k");
      }
    // If we have the "new" version, then parse it
    if( input.have_variable("Materials/"+material+"/ThermalConductivity/value") )
      {
        this->set_parameter
          (_k, input, "Materials/"+material+"/ThermalConductivity/value", _k);
      }
    // If instead we have the old version, use that.
    else if( input.have_variable("Materials/Conductivity/k") )
      {
        std::string warning = "WARNING: Specifying Materials/Conductivity/k is\n";
        warning += "         DEPRECATED. Please update to instead use\n";
        warning += "         Materials/MATERIAL_NAME/ThermalConductivity/value,\n";
        warning += "         where MATERIAL_NAME is given by Physics/PHYSICS_CLASS/material.\n";
        grins_warning(warning);

        this->set_parameter
          (_k, input, "Materials/Conductivity/k", _k);
      }

  }

  ConstantConductivity::~ConstantConductivity()
  {
    return;
  }

} // end namespace GRINS
