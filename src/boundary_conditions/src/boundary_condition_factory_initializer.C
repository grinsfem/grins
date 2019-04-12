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
#include "grins/boundary_condition_factory_initializer.h"

// GRINS
#include "grins/catalycity_factories.h"
#include "grins/constant_function_dirichlet_bc_factory.h"
#include "grins/gas_catalytic_wall_neumann_bc_factories.h"
#include "grins/homogeneous_dirichlet_bc_factory.h"
#include "grins/homogeneous_neumann_bc_factory.h"
#include "grins/parsed_function_dirichlet_bc_factory.h"
#include "grins/parsed_function_neumann_bc_factory.h"
#include "grins/symmetry_type_bc_factories.h"

// GRINS-Deprecated
#include "grins/catalycity_factories_old_style.h"
#include "grins/gas_catalytic_wall_neumann_bc_old_style_factories.h"
#include "grins/isothermal_dirichlet_old_style_bc_factory.h"
#include "grins/parsed_function_dirichlet_old_style_bc_factory.h"
#include "grins/parsed_function_neumann_old_style_bc_factory.h"
#include "grins/prescribed_vector_value_dirichlet_old_style_bc_factory.h"

namespace GRINS
{
  BoundaryConditionFactoryInitializer::BoundaryConditionFactoryInitializer()
  {
    static ConstantCatalycityFactory grins_factory_constant_catalycity("constant");
    static ArrheniusCatalycityFactory grins_factory_arrhenius_catalycity("arrhenius");
    static PowerLawCatalycityFactory grins_factory_power_law_catalycity("power");

    static ConstantFunctionDirichletBCFactory grins_factory_constant_dirichlet("constant_dirichlet");
    static ConstantFunctionDirichletBCFactory grins_factory_constant_displacement("constant_displacement");
    static ConstantFunctionDirichletBCFactory grins_factory_constant_isothermal("isothermal");
    static MoleFractionsDirichletBCFactory grins_factory_mole_fractions("mole_fractions");

    static GasRecombinationCatalyticWallNeumannBCFactory
      grins_factory_gas_recomb_catalytic_wall_neumann_bc("gas_recombination_catalytic_wall");

    static GasSolidCatalyticWallNeumannBCFactory
      grins_factory_gas_solid_catalytic_wall_neumann_bc("gas_solid_catalytic_wall");

    static HomogeneousDirichletBCFactory grins_factory_homogeneous_dirichlet("homogeneous_dirichlet");
    static HomogeneousDirichletBCFactory grins_factory_no_slip("no_slip");
    static HomogeneousDirichletBCFactory grins_factory_pinned("pinned");

    static HomogeneousNeumannBCFactory grins_factory_homogeneous_neumann("homogeneous_neumann");
    static HomogeneousNeumannBCFactory grins_factory_adiabatic("adiabatic");
    static HomogeneousNeumannBCFactory grins_factory_temp_axisym("Temperature_axisymmetric");
    static HomogeneousNeumannBCFactory grins_factory_species_axisym("SpeciesMassFractions_axisymmetric");

    static ParsedDirichletBCFactory grins_factory_parsed_dirichlet("parsed_dirichlet");
    static ParsedFEMDirichletBCFactory grins_factory_parsed_fem_dirichlet("parsed_fem_dirichlet");
    static ParsedDirichletBCFactory grins_factory_parsed_displacement("parsed_displacement");

    static ParsedFunctionNeumannBCFactory<libMesh::FunctionBase<libMesh::Number> >
      grins_factory_parsed_neumann("parsed_neumann");

    static ParsedFunctionNeumannBCFactory<libMesh::FEMFunctionBase<libMesh::Number> >
      grins_factory_parsed_fem_neumann("parsed_fem_neumann");

    static ParsedTractionBCFactory<libMesh::FunctionBase<libMesh::Number> >
      grins_factory_parsed_traction("parsed_traction");

    static ParsedTractionBCFactory<libMesh::FEMFunctionBase<libMesh::Number> >
      grins_factory_parsed_fem_traction("parsed_fem_traction");

    static ParsedTractionBCFactory<libMesh::FunctionBase<libMesh::Number> >
      grins_factory_constant_traction("constant_traction");

    static YZSymmetryBCFactory grins_factory_yz_symmetry("yz_symmetry");
    static XZSymmetryBCFactory grins_factory_xz_symmetry("xz_symmetry");
    static XYSymmetryBCFactory grins_factory_xy_symmetry("xy_symmetry");
    static RollerXBCFactory grins_factory_roller_x("roller_x");
    static RollerYBCFactory grins_factory_roller_y("roller_y");
    static RollerZBCFactory grins_factory_roller_z("roller_z");
    static AxisymmetryBCFactory grins_factory_velocity_axisymmetry("Velocity_axisymmetric");

    // Deprecated objects that will be removed in the future
    static ConstantCatalycityFactoryOldStyle grins_factory_constant_catalycity_old_style("constant_old_style");
    static ArrheniusCatalycityFactoryOldStyle grins_factory_arrhenius_catalycity_old_style("arrhenius_old_style");
    static PowerLawCatalycityFactoryOldStyle grins_factory_power_law_catalycity_old_style("power_old_style");

    static GasRecombinationCatalyticWallNeumannBCOldStyleFactory
      grins_factory_gas_recomb_catalytic_wall_neumann_bc_old_style("gas_recombination_catalytic_wall_old_style");

    static GasSolidCatalyticWallNeumannBCOldStyleFactory
      grins_factory_gas_solid_catalytic_wall_neumann_bc_old_style("gas_solid_catalytic_wall_old_style");

    static HomogeneousDirichletBCFactory grins_factory_no_slip_old_style("no_slip_old_style");
    static HomogeneousDirichletBCFactory grins_factory_pinned_old_style("pinned_old_style");

    static HomogeneousNeumannBCFactory grins_factory_adiabatic_old_style("adiabatic_old_style");
    static HomogeneousNeumannBCFactory grins_factory_adiabatic_wall_old_style("adiabatic_wall_old_style");
    static HomogeneousNeumannBCFactory grins_factory_temp_axisym_old_style("Temperature_axisymmetric_old_style");
    static HomogeneousNeumannBCFactory grins_factory_species_axisym_old_style("SpeciesMassFractions_axisymmetric_old_style");

    static IsothermalDirichletOldStyleBCFactory grins_factory_isothermal_wall_old_style("isothermal_wall_old_style");
    static IsothermalDirichletOldStyleBCFactory grins_factory_isothermal_old_style("isothermal_old_style");

    static ParsedDirichletOldStyleBCFactory grins_factory_parsed_dirichlet_old_style("parsed_dirichlet_old_style");
    static ParsedDirichletOldStyleBCFactory grins_factory_constant_dirichlet_old_style("constant_dirichlet_old_style");
    static ParsedFEMDirichletOldStyleBCFactory grins_factory_parsed_fem_dirichlet_old_style("parsed_fem_dirichlet_old_style");

    static TractionOldStyleBCFactory<libMesh::FunctionBase<libMesh::Number> >
      grins_factory_traction_old_style_functionbase("constant_traction_old_style");

    static PrescribedVelDirichletOldStyleBCFactory
      grins_factory_prescribed_vel_old_style("prescribed_vel_old_style");
    static PrescribedDispDirichletOldStyleBCFactory
      grins_factory_constant_displacement_old_style("constant_displacement_old_style");
    static PrescribedSpeciesDirichletOldStyleBCFactory
      grins_factory_prescribed_species_old_style("prescribed_species_old_style");
    static PrescribedMoleFractionsDirichletOldStyleBCFactory
      grins_factory_prescribed_mole_fraccs_old_style("prescribed_mole_fracs_old_style");

    static YZSymmetryBCFactory grins_factory_yz_symmetry_old_style("yz_symmetry_old_style");
    static XZSymmetryBCFactory grins_factory_xz_symmetry_old_style("xz_symmetry_old_style");
    static XYSymmetryBCFactory grins_factory_xy_symmetry_old_style("xy_symmetry_old_style");
    static RollerXBCFactory grins_factory_roller_x_old_style("roller_x_old_style");
    static RollerYBCFactory grins_factory_roller_y_old_style("roller_y_old_style");
    static RollerZBCFactory grins_factory_roller_z_old_style("roller_z_old_style");
    static AxisymmetryBCFactory grins_factory_velocity_axisymmetry_old_style("Velocity_axisymmetric_old_style");
  }
}// end namespace GRINS
