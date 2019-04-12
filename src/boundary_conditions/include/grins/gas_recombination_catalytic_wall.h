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

#ifndef GRINS_GAS_RECOMBINATION_CATALYTIC_WALL_H
#define GRINS_GAS_RECOMBINATION_CATALYTIC_WALL_H

// GRINS
#include "grins/var_typedefs.h"
#include "grins/catalytic_wall_base.h"

namespace GRINS
{

  template<typename Chemistry>
  class GasRecombinationCatalyticWall : public CatalyticWallBase<Chemistry>
  {
  public:

    GasRecombinationCatalyticWall( std::shared_ptr<Chemistry>& chem,
                                   std::shared_ptr<CatalycityBase>& gamma,
                                   const std::vector<VariableIndex>& species_vars,
                                   VariableIndex T_var,
                                   libMesh::Real p0,
                                   unsigned int reactant_species_idx,
                                   unsigned int product_species_idx );

    //! Deprecated
    GasRecombinationCatalyticWall( const Chemistry& chem_mixture,
                                   CatalycityBase& gamma,
                                   const unsigned int reactant_species_idx,
                                   const unsigned int product_species_idx );

    virtual ~GasRecombinationCatalyticWall(){};

    //! Deprecated
    virtual void init( const libMesh::FEMSystem& system );

    //! Deprecated
    virtual void apply_fluxes( AssemblyContext& context,
                               const CachedValues& cache,
                               const bool request_jacobian );

    virtual bool eval_flux( bool compute_jacobian,
                            AssemblyContext& context,
                            libMesh::Real sign,
                            bool is_axisymmetric );

    libMesh::Real compute_reactant_mass_flux( const libMesh::Real rho,
                                              const libMesh::Real Y_r,
                                              const libMesh::Real T );

    libMesh::Real compute_product_mass_flux( const libMesh::Real rho,
                                             const libMesh::Real Y_r,
                                             const libMesh::Real T );

  protected:

    const unsigned int _reactant_species_idx;
    VariableIndex _reactant_var_idx;

    const unsigned int _product_species_idx;
    VariableIndex _product_var_idx;

  };

  template<typename Chemistry>
  inline
  libMesh::Real GasRecombinationCatalyticWall<Chemistry>::compute_reactant_mass_flux( const libMesh::Real rho,
                                                                                      const libMesh::Real Y_r,
                                                                                      const libMesh::Real T )
  {
    const libMesh::Real rho_r = rho*Y_r;

    const libMesh::Real omega_dot = this->omega_dot( rho_r, T );

    return -omega_dot;
  }

  template<typename Chemistry>
  inline
  libMesh::Real GasRecombinationCatalyticWall<Chemistry>::compute_product_mass_flux( const libMesh::Real rho,
                                                                                     const libMesh::Real Y_r,
                                                                                     const libMesh::Real T )
  {
    return -(this->compute_reactant_mass_flux(rho,Y_r,T));
  }

} // end namespace GRINS

#endif // GRINS_GAS_RECOMBINATION_CATALYTIC_WALL_H
