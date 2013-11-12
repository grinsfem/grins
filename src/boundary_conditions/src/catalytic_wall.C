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


// This class
#include "grins/catalytic_wall.h"

// GRINS
#include "grins/cached_values.h"
#include "grins/constant_catalycity.h"
#include "grins/assembly_context.h"

namespace GRINS
{

  template<typename Chemistry>
  CatalyticWall<Chemistry>::CatalyticWall( const Chemistry& chemistry,
                                           const unsigned int species_index,
                                           CatalycityBase& gamma )
    : NeumannFuncObj(),
      _chemistry(chemistry),
      _species_index(species_index),
      _gamma_s( gamma.clone() ),
      _C( std::sqrt( chemistry.R(species_index)/(GRINS::Constants::two_pi) ) )
  {
    return;
  }

  template<typename Chemistry>
  CatalyticWall<Chemistry>::~CatalyticWall()
  {
    return;
  }

  template<typename Chemistry>
  void CatalyticWall<Chemistry>::init( const VariableIndex T_var )
  {
    _jac_vars.resize(1);
    _jac_vars[0] = T_var;

    return;
  }

  template<typename Chemistry>
  libMesh::Real CatalyticWall<Chemistry>::normal_value( const AssemblyContext& /*context*/,
                                                        const CachedValues& cache,
                                                        const unsigned int qp )
  {
    const libMesh::Real rho = cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp];
    
    const libMesh::Real w_s = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp][_species_index];
    
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];

    const libMesh::Real rho_s = rho*w_s;

    return this->omega_dot( rho_s, T );
  }

  template<typename Chemistry>
  libMesh::Real CatalyticWall<Chemistry>::normal_derivative( const AssemblyContext& /*context*/,
                                                             const CachedValues& cache,
                                                             const unsigned int qp )
  {
    const libMesh::Real rho = cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp];
    
    const libMesh::Real w_s = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp][_species_index];
    
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];

    const libMesh::Real rho_s = rho*w_s;

    const libMesh::Real R = _chemistry.R_mix( cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp] );
    
    return this->domega_dot_dws( rho_s, w_s, T, R );
  }

  template<typename Chemistry>
  libMesh::Real CatalyticWall<Chemistry>::normal_derivative( const AssemblyContext& /*context*/,
                                                             const CachedValues& cache,
                                                             const unsigned int qp, 
                                                             const GRINS::VariableIndex jac_var )
  {
    libmesh_assert_equal_to( jac_var, _jac_vars[0] );

    const libMesh::Real rho = cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp];
    
    const libMesh::Real w_s = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp][_species_index];
    
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];

    const libMesh::Real rho_s = rho*w_s;

    return this->domega_dot_dT( rho_s, T );
  }

} // end namespace GRINS
