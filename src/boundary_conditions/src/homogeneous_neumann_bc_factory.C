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
#include "grins/homogeneous_neumann_bc_factory.h"

namespace GRINS
{
  HomogeneousNeumannBCFactory grins_factory_homogeneous_neumann("homogeneous_neumann");
  HomogeneousNeumannBCFactory grins_factory_adiabatic("adiabatic");
  HomogeneousNeumannBCFactory grins_factory_adiabatic_old_style("adiabatic_old_style");
  HomogeneousNeumannBCFactory grins_factory_adiabatic_wall_old_style("adiabatic_wall_old_style");
  HomogeneousNeumannBCFactory grins_factory_temp_axisym("Temperature_axisymmetric");
  HomogeneousNeumannBCFactory grins_factory_temp_axisym_old_style("Temperature_axisymmetric_old_style");
  HomogeneousNeumannBCFactory grins_factory_species_axisym("SpeciesMassFractions_axisymmetric");
  HomogeneousNeumannBCFactory grins_factory_species_axisym_old_style("SpeciesMassFractions_axisymmetric_old_style");
} // end namespace GRINS
