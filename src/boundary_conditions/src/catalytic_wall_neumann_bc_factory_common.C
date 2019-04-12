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
#include "grins/catalytic_wall_neumann_bc_factory_common.h"

// GRINS
#include "grins/catalycity_factory_abstract.h"
#include "grins/fe_variables_base.h"
#include "grins/multicomponent_variable.h"
#include "grins/single_variable.h"
#include "grins/variable_warehouse.h"
#include "grins/gas_recombination_catalytic_wall_neumann_bc_factory_impl.h"
#include "grins/gas_solid_catalytic_wall_neumann_bc_factory_impl.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  template<typename ImpType>
  std::shared_ptr<NeumannBCAbstract>
  CatalyticWallNeumannBCFactoryCommon<ImpType>::build_catalytic_wall_common( const GetPot& input,
                                                                             const FEVariablesBase& fe_var,
                                                                             const std::string& material,
                                                                             const std::string& reaction,
                                                                             std::shared_ptr<CatalycityBase>& gamma_ptr,
                                                                             libMesh::Real p0,
                                                                             std::string& thermochem_lib )
  {
    // Extract species variables and material
    std::vector<VariableIndex> species_vars;
    this->extract_species_vars( fe_var, species_vars );

    libmesh_assert_equal_to( fe_var.active_var_names().size(), species_vars.size() );

    // Extract Temperature variable index.
    VariableIndex T_var = this->extract_temp_var();

    return _wall_impl.build_catalytic_wall(input,reaction,gamma_ptr,species_vars,
                                           material,T_var,p0,thermochem_lib);
  }

  template<typename ImpType>
  void CatalyticWallNeumannBCFactoryCommon<ImpType>::extract_species_vars
  ( const FEVariablesBase& fe_var, std::vector<VariableIndex>& species_vars ) const
  {
    species_vars = fe_var.var_indices();
  }

  template<typename ImpType>
  void CatalyticWallNeumannBCFactoryCommon<ImpType>::extract_material( const FEVariablesBase& fe_var,
                                                                       std::string& material ) const
  {
    const SpeciesMassFractionsVariable& species_fe_var =
      libMesh::cast_ref<const SpeciesMassFractionsVariable&>(fe_var);

    material = species_fe_var.material();
  }

  template<typename ImpType>
  VariableIndex CatalyticWallNeumannBCFactoryCommon<ImpType>::extract_temp_var() const
  {
    // This will throw an error is the temperature variables are not there
    const FEVariablesBase& temp_fe_var_base =
      GRINSPrivate::VariableWarehouse::get_variable(VariablesParsing::temperature_section());

    const PrimitiveTempFEVariables& temp_fe_var =
      libMesh::cast_ref<const PrimitiveTempFEVariables&>(temp_fe_var_base);

    std::vector<VariableIndex> temp_var = temp_fe_var.var_indices();

    libmesh_assert_equal_to( temp_var.size(), 1 );

    return temp_var[0];
  }

  template class CatalyticWallNeumannBCFactoryCommon<GasRecombinationCatalyticWallNeumannBCFactoryImpl>;
  template class CatalyticWallNeumannBCFactoryCommon<GasSolidCatalyticWallNeumannBCFactoryImpl>;

} // end namespace GRINS
