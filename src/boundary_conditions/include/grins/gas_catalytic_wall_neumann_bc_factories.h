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

#ifndef GRINS_GAS_CATALYTIC_WALL_NEUMANN_BC_FACTORIES_H
#define GRINS_GAS_CATALYTIC_WALL_NEUMANN_BC_FACTORIES_H

// GRINS
#include "grins/catalytic_wall_neumann_bc_factory_base.h"
#include "grins/gas_recombination_catalytic_wall_neumann_bc_factory_impl.h"
#include "grins/gas_solid_catalytic_wall_neumann_bc_factory_impl.h"

namespace GRINS
{
  class GasRecombinationCatalyticWallNeumannBCFactory :
    public CatalyticWallNeumannBCFactoryBase<GasRecombinationCatalyticWallNeumannBCFactoryImpl>
  {
  public:

    GasRecombinationCatalyticWallNeumannBCFactory( const std::string& bc_type_name )
      : CatalyticWallNeumannBCFactoryBase<GasRecombinationCatalyticWallNeumannBCFactoryImpl>(bc_type_name)
    {}

    ~GasRecombinationCatalyticWallNeumannBCFactory(){};

  };

  class GasSolidCatalyticWallNeumannBCFactory :
    public CatalyticWallNeumannBCFactoryBase<GasSolidCatalyticWallNeumannBCFactoryImpl>
  {
  public:

    GasSolidCatalyticWallNeumannBCFactory( const std::string& bc_type_name )
      : CatalyticWallNeumannBCFactoryBase<GasSolidCatalyticWallNeumannBCFactoryImpl>(bc_type_name)
    {}

    ~GasSolidCatalyticWallNeumannBCFactory(){};

  };
} // end namespace GRINS

#endif // GRINS_GAS_CATALYTIC_WALL_NEUMANN_BC_FACTORIES_H
