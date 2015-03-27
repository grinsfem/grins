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
#include "grins/grins_physics_names.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  void PhysicsFactoryHelper::parse_viscosity_model( const GetPot& input,
                                                    const std::string& physics,
                                                    std::string& model )
  {
    model = input( "Physics/"+incompressible_navier_stokes+"/viscosity_model", "constant" );
    return;
  }

  void PhysicsFactoryHelper::parse_conductivity_model( const GetPot& input,
                                                       const std::string& physics,
                                                       std::string& model )
  {
    model = input( "Physics/"+heat_transfer+"/conductivity_model", "constant" );
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
