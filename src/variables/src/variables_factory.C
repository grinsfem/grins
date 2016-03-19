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
#include "grins/variables_factory.h"

// GRINS
#include "grins/variables_parsing.h"
#include "grins/displacement_variables.h"
#include "grins/pressure_variable.h"
#include "grins/primitive_temp_variables.h"
#include "grins/thermo_pressure_variable.h"
#include "grins/turbulence_variables.h"
#include "grins/velocity_variables.h"
#include "grins/species_mass_fracs_variables.h"

namespace GRINS
{
  template<typename DerivedVariables>
  libMesh::UniquePtr<VariablesBase> VariablesFactory<DerivedVariables>::build_vars( const GetPot& input )
  {
    return libMesh::UniquePtr<VariablesBase>( new DerivedVariables(input) );
  }

  // Instantiate all the "Basic" Variables factories.
  // These shouldn't be directly used by the user, we just need to instantiate them.
  VariablesFactory<DisplacementVariables> grins_factory_displacement_vars(VariablesParsing::displacement_section());
  VariablesFactory<PressureVariable> grins_factory_pressure_var(VariablesParsing::pressure_section());
  VariablesFactory<PrimitiveTempVariables> grins_factory_prim_temp_vars(VariablesParsing::temperature_section());
  VariablesFactory<ThermoPressureVariable> grins_factory_thermo_pressure_var(VariablesParsing::thermo_pressure_section());
  VariablesFactory<VelocityVariables> grins_factory_velocity_vars(VariablesParsing::velocity_section());
  VariablesFactory<SpeciesMassFractionsVariables> grins_factory_species_mass_frac_vars(VariablesParsing::species_mass_fractions_section());
} // end namespace GRINS
