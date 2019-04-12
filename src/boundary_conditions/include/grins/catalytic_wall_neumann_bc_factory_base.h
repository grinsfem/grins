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

#ifndef GRINS_CATALYTIC_WALL_NEUMANN_BC_FACTORY_BASE_H
#define GRINS_CATALYTIC_WALL_NEUMANN_BC_FACTORY_BASE_H

// GRINS
#include "grins/neumann_bc_factory_abstract.h"
#include "grins/catalytic_wall_neumann_bc_factory_common.h"

namespace GRINS
{

  template<typename ImplType>
  class CatalyticWallNeumannBCFactoryBase : public NeumannBCFactoryAbstract,
                                            public CatalyticWallNeumannBCFactoryCommon<ImplType>
  {
  public:

    CatalyticWallNeumannBCFactoryBase( const std::string& bc_type_name )
      : NeumannBCFactoryAbstract(bc_type_name),
        CatalyticWallNeumannBCFactoryCommon<ImplType>()
    {}

    ~CatalyticWallNeumannBCFactoryBase(){};

  protected:

    virtual std::shared_ptr<NeumannBCAbstract>
    build_neumann_func( const GetPot& input,
                        MultiphysicsSystem& system,
                        const FEVariablesBase& fe_var,
                        const std::string& section );

    //! Parse the reaction.
    std::string parse_reaction( const GetPot& input, const std::string& section ) const;

    std::string parse_thermochem_model( const GetPot& input,
                                        const std::string& material ) const;

    libMesh::Real parse_thermo_pressure( const GetPot& input,
                                         const std::string& material ) const;

    std::shared_ptr<CatalycityBase>
    build_catalycity( const GetPot& input, const std::string& section ) const;
  };
} // end namespace GRINS

#endif // GRINS_CATALYTIC_WALL_NEUMANN_BC_FACTORY_BASE_H
