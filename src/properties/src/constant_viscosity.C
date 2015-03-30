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
#include "grins/constant_viscosity.h"

//GRINS
#include "grins/grins_physics_names.h"
#include "grins/common.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{

  ConstantViscosity::ConstantViscosity( const GetPot& input )
    : ParameterUser("ConstantViscosity"),
      _mu(1.0)
  {
    // Warning about this constructor being deprecated
    {
      std::string warning = "WARNING: Use of this constructor is DEPRECATED.\n";
      warning += "         Please update to use constructor with input material name.\n";
      grins_warning(warning);
    }

    if( !input.have_variable("Materials/Viscosity/mu") )
      {
        libmesh_warning("No Materials/Viscosity/mu specified!\n");

	// Try and get the viscosity from other specifications
        this->set_parameter
	  (_mu, input,
           "Physics/"+incompressible_navier_stokes+"/mu", _mu);
	
      }
    else
      this->set_parameter
        (_mu, input, "Materials/Viscosity/mu", _mu);
  }

  ConstantViscosity::~ConstantViscosity()
  {
    return;
  }

} // namespace GRINS
