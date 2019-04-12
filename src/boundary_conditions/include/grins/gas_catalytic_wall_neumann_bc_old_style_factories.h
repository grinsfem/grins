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

#ifndef GRINS_GAS_CATALYTIC_WALL_NEUMANN_BC_OLD_STYLE_FACTORIES_H
#define GRINS_GAS_CATALYTIC_WALL_NEUMANN_BC_OLD_STYLE_FACTORIES_H

// GRINS
#include "grins/catalytic_wall_neumann_bc_old_style_factory_base.h"
#include "grins/gas_recombination_catalytic_wall_neumann_bc_factory_impl.h"
#include "grins/gas_solid_catalytic_wall_neumann_bc_factory_impl.h"

namespace GRINS
{
  class GasRecombinationCatalyticWallNeumannBCOldStyleFactory :
    public CatalyticWallNeumannBCOldStyleFactoryBase<GasRecombinationCatalyticWallNeumannBCFactoryImpl>
  {
  public:

    GasRecombinationCatalyticWallNeumannBCOldStyleFactory( const std::string& bc_type_name )
      : CatalyticWallNeumannBCOldStyleFactoryBase<GasRecombinationCatalyticWallNeumannBCFactoryImpl>(bc_type_name)
    {}

    ~GasRecombinationCatalyticWallNeumannBCOldStyleFactory(){};

  protected:

    virtual std::string reactant_for_catalycity(const std::string& reaction) const
    {
      std::string reactant;
      std::string product; // dummy in this function
      this->_wall_impl.parse_reactant_and_product(reaction,reactant,product);
      return reactant;
    }

    virtual std::string catalytic_wall_prefix_str() const
    {
      return "wall_catalytic_reactions";
    }
  };

  class GasSolidCatalyticWallNeumannBCOldStyleFactory :
    public CatalyticWallNeumannBCOldStyleFactoryBase<GasSolidCatalyticWallNeumannBCFactoryImpl>
  {
  public:

    GasSolidCatalyticWallNeumannBCOldStyleFactory( const std::string& bc_type_name )
      : CatalyticWallNeumannBCOldStyleFactoryBase<GasSolidCatalyticWallNeumannBCFactoryImpl>(bc_type_name)
    {}

    ~GasSolidCatalyticWallNeumannBCOldStyleFactory(){};

    virtual std::string reactant_for_catalycity(const std::string& reaction) const
    {
      std::string gas_reactant;
      std::string solid_reactant; // dummy in this function
      std::string product; // dummy in this function
      this->_wall_impl.parse_reactants_and_product(reaction,gas_reactant,solid_reactant,product);
      return gas_reactant;
    }

    virtual std::string catalytic_wall_prefix_str() const
    {
      return "wall_gas_solid_reactions";
    }

  };
} // end namespace GRINS

#endif // GRINS_GAS_CATALYTIC_WALL_NEUMANN_BC_OLD_STYLE_FACTORIES_H
