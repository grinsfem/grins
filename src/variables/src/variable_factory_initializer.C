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
#include "grins/variable_factory_initializer.h"

// GRINS
#include "grins/variable_factory.h"
#include "grins/single_variable.h"
#include "grins/multicomponent_variable.h"
#include "grins/multi_component_vector_variable.h"

namespace GRINS
{
  VariableFactoryInitializer::VariableFactoryInitializer()
  {
    static VariableFactoryBasic<DisplacementVariable>
      grins_factory_disp_fe_var(VariablesParsing::displacement_section());

    static VariableFactoryBasic<SingleVariable>
      grins_factory_single_var(VariablesParsing::single_var_section());

    static VariableFactoryBasic<PressureFEVariable>
      grins_factory_press_fe_var(VariablesParsing::pressure_section());

    static VariableFactoryBasic<PrimitiveTempFEVariables>
      grins_factory_temp_fe_var(VariablesParsing::temperature_section());

    static SpeciesVariableFactory<SpeciesMassFractionsVariable>
      grins_factory_species_mass_frac_fe_var(VariablesParsing::species_mass_fractions_section());

    static ScalarVariableFactory<ThermoPressureVariable>
      grins_factory_thermo_press_fe_var(VariablesParsing::thermo_pressure_section());

    static VariableFactoryBasic<TurbulenceFEVariables>
      grins_factory_turb_fe_var(VariablesParsing::turbulence_section());

    static VariableFactoryBasic<VelocityVariable>
      grins_factory_velocity_fe_var(VariablesParsing::velocity_section());

    static ScalarVariableFactory<ScalarVariable>
      grins_factory_scalar_var(VariablesParsing::scalar_var_section());
  }
} // end namespace GRINS
