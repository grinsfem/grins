//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "grins/catalytic_wall.h"

// GRINS
#include "grins/cached_values.h"
#include "grins/chemical_mixture.h"

// libMesh
#include "libmesh/fem_context.h"

namespace GRINS
{

  CatalyticWall::CatalyticWall( const ChemicalMixture& chem_mixture,
				const unsigned int species_index,
				const VariableIndex T_var,
				const libMesh::Real gamma )
    : NeumannFuncObj(),
      _chem_mixture(chem_mixture),
      _species_index(species_index),
      _T_var(T_var),
      _helper( _chem_mixture.R(_species_index), _chem_mixture.M(_species_index), gamma )
  {
    _jac_vars.resize(1);
    _jac_vars[0] = _T_var;

    return;
  }

  CatalyticWall::~CatalyticWall()
  {
    return;
  }

  libMesh::Real CatalyticWall::normal_value( const libMesh::FEMContext& /*context*/,
					     const CachedValues& cache,
					     const unsigned int qp )
  {
    const libMesh::Real rho = cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp];
    
    const libMesh::Real w_s = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp][_species_index];
    
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];

    const libMesh::Real rho_s = rho*w_s;

    return this->omega_dot( rho_s, T );
  }

  libMesh::Real CatalyticWall::normal_derivative( const libMesh::FEMContext& /*context*/,
						  const CachedValues& cache,
						  const unsigned int qp )
  {
    const libMesh::Real rho = cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp];
    
    const libMesh::Real w_s = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp][_species_index];
    
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];

    const libMesh::Real rho_s = rho*w_s;

    const libMesh::Real R = _chem_mixture.R( cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp] );
    
    return this->domega_dot_dws( rho_s, w_s, T, R );
  }

  libMesh::Real CatalyticWall::normal_derivative( const libMesh::FEMContext& /*context*/,
						  const CachedValues& cache,
						  const unsigned int qp, 
						  const GRINS::VariableIndex jac_var )
  {
    libmesh_assert_equal_to( jac_var, _T_var );

    const libMesh::Real rho = cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp];
    
    const libMesh::Real w_s = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp][_species_index];
    
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];

    const libMesh::Real rho_s = rho*w_s;

    return this->domega_dot_dT( rho_s, T );
  }

} // namespace GRINS
