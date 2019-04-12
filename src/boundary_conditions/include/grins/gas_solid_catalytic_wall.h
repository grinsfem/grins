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

#ifndef GRINS_GAS_SOLID_CATALYTIC_WALL_H
#define GRINS_GAS_SOLID_CATALYTIC_WALL_H

// GRINS
#include "grins/var_typedefs.h"
#include "grins/catalytic_wall_base.h"

namespace GRINS
{

  template<typename Chemistry>
  class GasSolidCatalyticWall : public CatalyticWallBase<Chemistry>
  {
  public:

    GasSolidCatalyticWall(std::shared_ptr<Chemistry>& chem,
                          std::shared_ptr<CatalycityBase>& gamma,
                          const std::vector<VariableIndex>& species_vars,
                          VariableIndex T_var,
                          libMesh::Real p0,
                          unsigned int reactant_gas_species_idx,
                          unsigned int reactant_solid_species_idx,
                          unsigned int product_species_idx);

    //! Deprecated
    GasSolidCatalyticWall( const Chemistry& chem_mixture,
                           CatalycityBase& gamma,
                           const unsigned int reactant_gas_species_idx,
                           const unsigned int reactant_solid_species_idx,
                           const unsigned int product_species_idx );

    virtual ~GasSolidCatalyticWall(){};

    //! Deprecated
    virtual void init( const libMesh::FEMSystem& system );

    virtual bool eval_flux( bool compute_jacobian,
                            AssemblyContext& context,
                            libMesh::Real sign,
                            bool is_axisymmetric );

    //! Deprecated
    virtual void apply_fluxes( AssemblyContext& context,
                               const CachedValues& cache,
                               const bool request_jacobian );

    libMesh::Real compute_reactant_gas_mass_flux( const libMesh::Real rho,
                                                  const libMesh::Real Y_r,
                                                  const libMesh::Real T );

    // kg/m^2-s
    libMesh::Real compute_reactant_solid_mass_consumption( const libMesh::Real rho,
                                                           const libMesh::Real Y_r,
                                                           const libMesh::Real T );

    libMesh::Real compute_product_mass_flux( const libMesh::Real rho,
                                             const libMesh::Real Y_r,
                                             const libMesh::Real T );

    libMesh::Real compute_reactant_solid_mass_consumption_dT( const libMesh::Real rho,
                                                              const libMesh::Real Y_r,
                                                              const libMesh::Real T );

    libMesh::Real compute_reactant_solid_mass_consumption_dYs( const libMesh::Real rho,
                                                               const std::vector<libMesh::Real> Y,
                                                               const libMesh::Real T );

  protected:

    const unsigned int _reactant_gas_species_idx;
    VariableIndex _reactant_gas_var_idx;

    const unsigned int _reactant_solid_species_idx;

    const unsigned int _product_species_idx;
    VariableIndex _product_var_idx;

  };

  template<typename Chemistry>
  inline
  libMesh::Real GasSolidCatalyticWall<Chemistry>::compute_reactant_gas_mass_flux( const libMesh::Real rho,
                                                                                  const libMesh::Real Y_r,
                                                                                  const libMesh::Real T )
  {
    const libMesh::Real rho_r = rho*Y_r;

    const libMesh::Real omega_dot = this->omega_dot( rho_r, T );

    return -omega_dot;
  }

  template<typename Chemistry>
  inline
  libMesh::Real GasSolidCatalyticWall<Chemistry>::compute_reactant_solid_mass_consumption( const libMesh::Real rho,
                                                                                           const libMesh::Real Y_r,
                                                                                           const libMesh::Real T )
  {
    const libMesh::Real rho_r = rho*Y_r;

    const libMesh::Real M_r = this->_chemistry.M(_reactant_gas_species_idx);

    const libMesh::Real M_solid = this->_chemistry.M(_reactant_solid_species_idx);

    const libMesh::Real omega_dot = this->omega_dot( rho_r, T )*M_solid/M_r;

    return -omega_dot;
  }

  template<typename Chemistry>
  inline
  libMesh::Real GasSolidCatalyticWall<Chemistry>::compute_product_mass_flux( const libMesh::Real rho,
                                                                             const libMesh::Real Y_r,
                                                                             const libMesh::Real T )
  {
    const libMesh::Real rho_r = rho*Y_r;

    const libMesh::Real M_r = this->_chemistry.M(_reactant_gas_species_idx);

    const libMesh::Real M_p = this->_chemistry.M(_product_species_idx);

    const libMesh::Real omega_dot = this->omega_dot( rho_r, T )*M_p/M_r;

    return omega_dot;
  }

  template<typename Chemistry>
  inline
  libMesh::Real GasSolidCatalyticWall<Chemistry>::compute_reactant_solid_mass_consumption_dT( const libMesh::Real rho,
                                                                                              const libMesh::Real Y_r,
                                                                                              const libMesh::Real T )
  {
    const libMesh::Real rho_r = rho*Y_r;

    const libMesh::Real M_r = this->_chemistry.M(_reactant_gas_species_idx);

    const libMesh::Real M_solid = this->_chemistry.M(_reactant_solid_species_idx);

    const libMesh::Real domega_dot_dT = this->domega_dot_dT( rho_r, T )*M_solid/M_r;

    return -domega_dot_dT;
  }

  template<typename Chemistry>
  inline
  libMesh::Real GasSolidCatalyticWall<Chemistry>::compute_reactant_solid_mass_consumption_dYs( const libMesh::Real rho,
                                                                                               const std::vector<libMesh::Real> Y,
                                                                                               const libMesh::Real T )
  {
    const libMesh::Real Y_r = Y[_reactant_gas_species_idx];

    const libMesh::Real rho_r = rho*Y_r;

    const libMesh::Real M_r = this->_chemistry.M(_reactant_gas_species_idx);

    const libMesh::Real M_solid = this->_chemistry.M(_reactant_solid_species_idx);

    const libMesh::Real R = this->_chemistry.R_mix(Y);

    const libMesh::Real domega_dot_dYs = this->domega_dot_dws( rho_r, Y_r, T, R )*M_solid/M_r;

    return -domega_dot_dYs;
  }

} // end namespace GRINS

#endif // GRINS_GAS_SOLID_CATALYTIC_WALL_H
