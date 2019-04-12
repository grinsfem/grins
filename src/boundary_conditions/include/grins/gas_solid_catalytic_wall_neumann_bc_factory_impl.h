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

#ifndef GRINS_GAS_SOLID_CATALYTIC_WALL_NEUMANN_BC_IMPL_H
#define GRINS_GAS_SOLID_CATALYTIC_WALL_NEUMANN_BC_IMPL_H

// GRINS
#include "grins/gas_solid_catalytic_wall.h"
#include "grins/fe_variables_base.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  class GasSolidCatalyticWallNeumannBCFactoryImpl
  {
  public:

    GasSolidCatalyticWallNeumannBCFactoryImpl(){}

    ~GasSolidCatalyticWallNeumannBCFactoryImpl(){}

    std::shared_ptr<NeumannBCAbstract>
    build_catalytic_wall( const GetPot& input,
                          const std::string& reaction,
                          std::shared_ptr<CatalycityBase>& gamma_ptr,
                          const std::vector<VariableIndex>& species_vars,
                          const std::string& material,
                          VariableIndex T_var,
                          libMesh::Real p0,
                          const std::string& thermochem_lib );

    void parse_reactants_and_product( const std::string& reaction,
                                      std::string& gas_reactant,
                                      std::string& solid_reactant,
                                      std::string& product ) const;

  protected:

    template<typename ChemistryType>
    void build_wall_ptr( std::shared_ptr<ChemistryType> & chem_ptr,
                         std::shared_ptr<CatalycityBase> & gamma_ptr,
                         const std::string & gas_reactant,
                         const std::string & solid_reactant,
                         const std::string & product,
                         const std::vector<VariableIndex> & species_vars,
                         VariableIndex T_var,
                         libMesh::Real p0,
                         std::shared_ptr<NeumannBCAbstract> & catalytic_wall )
    {
      catalytic_wall.reset( new GasSolidCatalyticWall<ChemistryType>
                            ( chem_ptr,
                              gamma_ptr,
                              species_vars,
                              T_var,
                              p0,
                              chem_ptr->species_index(gas_reactant),
                              chem_ptr->species_index(solid_reactant),
                              chem_ptr->species_index(product) ) );
    }

  };

} // end namespace GRINS

#endif // GRINS_GAS_SOLID_CATALYTIC_WALL_NEUMANN_BC_IMPL_H
