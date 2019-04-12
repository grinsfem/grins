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

#ifndef GRINS_CATALYTIC_WALL_NEUMANN_BC_FACTORY_COMMON_H
#define GRINS_CATALYTIC_WALL_NEUMANN_BC_FACTORY_COMMON_H

// GRINS
#include "grins/catalycity_base.h"
#include "grins/var_typedefs.h"
#include "grins/neumann_bc_abstract.h"

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  // Forward declarations
  class FEVariablesBase;

  //! Factory helper class for building catalytic wall Neumann boundary conditions
  /*! Note for catalytic walls, we're currently assuming that both SpeciesMassFractions
    and Temperature FEVariables are in the system. */
  template<typename ImplType>
  class CatalyticWallNeumannBCFactoryCommon
  {
  public:

    CatalyticWallNeumannBCFactoryCommon(){};

    ~CatalyticWallNeumannBCFactoryCommon(){};

  protected:

    std::shared_ptr<NeumannBCAbstract>
    build_catalytic_wall_common( const GetPot& input,
                                 const FEVariablesBase& fe_var,
                                 const std::string& material,
                                 const std::string& reaction,
                                 std::shared_ptr<CatalycityBase>& gamma_ptr,
                                 libMesh::Real p0,
                                 std::string& thermochem_lib );

    ImplType _wall_impl;

    void extract_species_vars( const FEVariablesBase& fe_var,
                               std::vector<VariableIndex>& species_vars ) const;

    void extract_material( const FEVariablesBase& fe_var, std::string& material ) const;

    VariableIndex extract_temp_var() const;

  };

} // end namespace GRINS

#endif // GRINS_CATALYTIC_WALL_NEUMANN_BC_FACTORY_COMMON_H
