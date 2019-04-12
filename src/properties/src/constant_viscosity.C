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
#include "grins/constant_viscosity.h"

//GRINS
#include "grins/physics_naming.h"
#include "grins/common.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{

  ConstantViscosity::ConstantViscosity( const GetPot& input )
    : ParameterUser("ConstantViscosity"),
      ViscosityBase(),
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
          (_mu, input, "Physics/"+PhysicsNaming::incompressible_navier_stokes()+"/mu", _mu);
      }
    else
      {
        this->set_parameter
          (_mu, input, "Materials/Viscosity/mu", _mu);
      }
  }

  ConstantViscosity::ConstantViscosity( const GetPot& input, const std::string& material )
    : ParameterUser("ConstantViscosity"),
      ViscosityBase(),
      _mu( 0.0 ) // Initialize to nonsense value
  {
    this->check_input_consistency(input,material);

    if( input.have_variable("Materials/"+material+"/Viscosity/value") &&
        input.have_variable("Physics/"+PhysicsNaming::incompressible_navier_stokes()+"/mu") )
      {
        libmesh_error_msg("Error: Cannot specify both Materials/"+material+"/Viscosity/value and Physics/"+PhysicsNaming::incompressible_navier_stokes()+"/mu");
      }

    // If we have the "new" version, then parse it
    if( input.have_variable("Materials/"+material+"/Viscosity/value") )
      {
        this->set_parameter
          (_mu, input, "Materials/"+material+"/Viscosity/value", _mu);
      }
    // If instead we have the old version, use that.
    else if( input.have_variable("Materials/Viscosity/mu") )
      {
        this->old_mu_warning();

        this->set_parameter
          (_mu, input, "Materials/Viscosity/mu", _mu);
      }
    /* If we don't have the new version of materials parsing or
       explicitly have the old version, we're assuming an even older
       version. Both of the older versions are deprecated. */
    else if( !input.have_variable("Materials/Viscosity/mu") &&
             !input.have_variable("Materials/"+material+"/Viscosity/value") )
      {
        std::string warning = "WARNING: No Materials/Viscosity/mu or\n";
        warning += "       Materials/"+material+"/Viscosity/value specified!\n";
        warning += "       We are assuming then that you want to specify through\n";
        warning += "       Physics/"+PhysicsNaming::incompressible_navier_stokes()+"/mu.\n";
        warning += "       This is DEPRECATED. Please updated to use Materials/"+material+"/Viscosity/value.\n";
        grins_warning(warning);

        // Try and get the viscosity from other specifications
        this->set_parameter
          (_mu, input, "Physics/"+PhysicsNaming::incompressible_navier_stokes()+"/mu", _mu);
      }
    else
      {
        // This shouldn't happen
        libmesh_error();
      }

    // We'd better have postive viscosity when we're all done.
    libmesh_assert_greater( _mu, 0.0 );
  }

  ConstantViscosity::~ConstantViscosity()
  {
    return;
  }

} // namespace GRINS
