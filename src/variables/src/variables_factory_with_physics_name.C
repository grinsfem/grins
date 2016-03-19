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
#include "grins/variables_factory_with_physics_name.h"

// GRINS
#include "grins/variables_parsing.h"
#include "grins/generic_variable.h"

namespace GRINS
{
  template<typename DerivedVariables>
  libMesh::UniquePtr<VariablesBase> VariablesFactoryWithPhysicsName<DerivedVariables>::build_vars( const GetPot& input )
  {
    // Make sure user set the physics name
    if( _physics_name == std::string("DIE!") )
      libmesh_error_msg("ERROR: must call set_physics_name() before building VariablesBase!");

    libMesh::UniquePtr<VariablesBase> new_var( new DerivedVariables(input,_physics_name) );

    // Reset the _physics_name for error checking
    _physics_name = std::string("DIE!");

    return new_var;
  }

  // Instantiate all the Variables factories that require a physics_name in the constructor
  // These shouldn't be directly used by the user, we just need to instantiate them.
  VariablesFactoryWithPhysicsName<GenericVariable> grins_factory_generic_var(VariablesParsing::generic_section());

} // end namespace GRINS
