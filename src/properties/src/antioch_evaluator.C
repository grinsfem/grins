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


#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// C++
#include <limits>

// This class
#include "grins/antioch_evaluator.h"

// GRINS
#include "grins/antioch_mixture.h"
#include "grins/cached_values.h"

namespace GRINS
{
  template<typename Thermo>
  AntiochEvaluator<Thermo>::AntiochEvaluator( const AntiochMixture& mixture )
    : _chem( mixture ),
      _thermo( NULL ),
      _kinetics( new AntiochKinetics(mixture) ),
      _temp_cache( new Antioch::TempCache<libMesh::Real>(1.0) )
  {
    this->build_thermo( mixture );
    return;
  }

  template<typename Thermo>
  AntiochEvaluator<Thermo>::~AntiochEvaluator()
  {
    return;
  }

  template<typename Thermo>
  void AntiochEvaluator<Thermo>::omega_dot( const CachedValues& cache, unsigned int qp,
                                            std::vector<libMesh::Real>& omega_dot )
  {
    const libMesh::Real& T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const libMesh::Real rho = cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp];
    const libMesh::Real R_mix = cache.get_cached_values(Cache::MIXTURE_GAS_CONSTANT)[qp];
    const std::vector<libMesh::Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];

    this->check_and_reset_temp_cache(T);

    _kinetics->omega_dot( *(_temp_cache.get()), rho, R_mix, Y, omega_dot );

    return;
  }

  template<typename Thermo>
  void AntiochEvaluator<Thermo>::check_and_reset_temp_cache( const libMesh::Real& T )
  {
    if( _temp_cache->T != T )
      {
        _temp_cache.reset( new Antioch::TempCache<libMesh::Real>(T) );
      }

    return;
  }

  template<>
  libMesh::Real AntiochEvaluator<Antioch::CEAEvaluator<libMesh::Real> >::cp( const CachedValues& cache,
                                                                             unsigned int qp )
  {
    const libMesh::Real& T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const std::vector<libMesh::Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];

    this->check_and_reset_temp_cache(T);

    return _thermo->cp( *(_temp_cache.get()), Y );
  }

  template<>
  libMesh::Real AntiochEvaluator<Antioch::StatMechThermodynamics<libMesh::Real> >::cp( const CachedValues& cache,
                                                                                       unsigned int qp )
  {
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const std::vector<libMesh::Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];

    return _thermo->cp( T, T, Y );
  }

  template<>
  libMesh::Real AntiochEvaluator<Antioch::CEAEvaluator<libMesh::Real> >::cp( const libMesh::Real& T,
                                                                             const std::vector<libMesh::Real>& Y )
  {
    this->check_and_reset_temp_cache(T);

    return _thermo->cp( *(_temp_cache.get()), Y );
  }

  template<>
  libMesh::Real AntiochEvaluator<Antioch::StatMechThermodynamics<libMesh::Real> >::cp( const libMesh::Real& T,
                                                                                       const std::vector<libMesh::Real>& Y )
  {
    return _thermo->cp( T, T, Y );
  }

  template<>
  libMesh::Real AntiochEvaluator<Antioch::CEAEvaluator<libMesh::Real> >::cv( const CachedValues& cache,
                                                                             unsigned int qp )
  {
    const libMesh::Real& T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const std::vector<libMesh::Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];

    this->check_and_reset_temp_cache(T);

    return _thermo->cv( *(_temp_cache.get()), Y );
  }

  template<>
  libMesh::Real AntiochEvaluator<Antioch::StatMechThermodynamics<libMesh::Real> >::cv( const CachedValues& cache,
                                                                                       unsigned int qp )
  {
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const std::vector<libMesh::Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];

    return _thermo->cv( T, T, Y );
  }
     
  template<>
  libMesh::Real AntiochEvaluator<Antioch::CEAEvaluator<libMesh::Real> >::h_s( const CachedValues& cache,
                                                                              unsigned int qp,
                                                                              unsigned int species )
  {
    const libMesh::Real& T = cache.get_cached_values(Cache::TEMPERATURE)[qp];

    this->check_and_reset_temp_cache(T);

    return _thermo->h( *(_temp_cache.get()), species );
  }

  template<>
  libMesh::Real AntiochEvaluator<Antioch::StatMechThermodynamics<libMesh::Real> >::h_s( const CachedValues& cache,
                                                                                        unsigned int qp,
                                                                                        unsigned int species )
  {
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];

    return _thermo->h_tot( species, T ) + _chem.h_stat_mech_ref_correction(species);
  }

  template<>
  void AntiochEvaluator<Antioch::CEAEvaluator<libMesh::Real> >::h_s( const CachedValues& cache,
                                                                     unsigned int qp,
                                                                     std::vector<libMesh::Real>& h_s )
  {
    const libMesh::Real& T = cache.get_cached_values(Cache::TEMPERATURE)[qp];

    this->check_and_reset_temp_cache(T);

    _thermo->h( *(_temp_cache.get()), h_s );

    return;
  }

  template<>
  void AntiochEvaluator<Antioch::StatMechThermodynamics<libMesh::Real> >::h_s( const CachedValues& cache,
                                                                               unsigned int qp,
                                                                               std::vector<libMesh::Real>& h_s )
  {
    for( unsigned int s = 0; s < _chem.n_species(); s++ )
      {
        h_s[s] = this->h_s(cache,qp,s);
      }

    return;
  }

} // end namespace GRINS

#endif //GRINS_HAVE_ANTIOCH
