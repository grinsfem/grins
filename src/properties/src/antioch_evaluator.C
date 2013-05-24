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

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// This class
#include "grins/antioch_evaluator.h"

// GRINS
#include "grins/antioch_mixture.h"
#include "grins/cached_values.h"

namespace GRINS
{
  template<typename AntiochThermo>
  AntiochEvaluator<AntiochThermo>::AntiochEvaluator( const AntiochMixture& mixture )
    : _chem( mixture ),
      _thermo( mixture ),
      _kinetics( AntiochKinetics(mixture) ),
      _temp_cache( new Antioch::TempCache<libMesh::Real>(0.0) )
  {
    return;
  }

  template<typename AntiochThermo>
  AntiochEvaluator<AntiochThermo>::~AntiochEvaluator()
  {
    return;
  }

  template<typename AntiochThermo>
  libMesh::Real AntiochEvaluator<AntiochThermo>::cp( const CachedValues& cache, unsigned int qp )
  {
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const std::vector<libMesh::Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];

    /*! \todo AntiochStatMechThermo doesn't need the full TempCache, just T.
      Should we try and optimize for that case? */
    this->check_and_reset_temp_cache(T);

    return _thermo.cp( *(_temp_cache.get()), Y );
  }

  template<typename AntiochThermo>
  libMesh::Real AntiochEvaluator<AntiochThermo>::cv( const CachedValues& cache, unsigned int qp )
  {
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const std::vector<libMesh::Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];

    /*! \todo AntiochStatMechThermo doesn't need the full TempCache, just T.
      Should we try and optimize for that case? */
    this->check_and_reset_temp_cache(T);

    return _thermo.cv( *(_temp_cache.get()), Y );
  }
     
  template<typename AntiochThermo>
  libMesh::Real AntiochEvaluator<AntiochThermo>::h_s( const CachedValues& cache, unsigned int qp,
                                                      unsigned int species )
  {
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];

    /*! \todo AntiochStatMechThermo doesn't need the full TempCache, just T.
      Should we try and optimize for that case? */
    this->check_and_reset_temp_cache(T);

    return _thermo.h_s( *(_temp_cache.get()), species );
  }
  
  template<typename AntiochThermo>
  void AntiochEvaluator<AntiochThermo>::h_s( const CachedValues& cache, unsigned int qp,
                                             std::vector<libMesh::Real>& h_s )
  {
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];

    /*! \todo AntiochStatMechThermo doesn't need the full TempCache, just T.
      Should we try and optimize for that case? */
    this->check_and_reset_temp_cache(T);

    _thermo.h_s( *(_temp_cache.get()), h_s );
    
    return;
  }

  template<typename AntiochThermo>
  libMesh::Real AntiochEvaluator<AntiochThermo>::mu( const CachedValues& cache, unsigned int qp )
  {
    libmesh_not_implemented();
    return 0.0;
  }

  template<typename AntiochThermo>
  libMesh::Real AntiochEvaluator<AntiochThermo>::k( const CachedValues& cache, unsigned int qp )
  {
    libmesh_not_implemented();
    return 0.0;
  }

  template<typename AntiochThermo>
  void AntiochEvaluator<AntiochThermo>::D( const CachedValues& cache, unsigned int qp,
                                           std::vector<libMesh::Real>& D )
  {
    libmesh_not_implemented();
    return;
  }

  template<typename AntiochThermo>
  void AntiochEvaluator<AntiochThermo>::omega_dot( const CachedValues& cache, unsigned int qp,
                                                   std::vector<libMesh::Real>& omega_dot )
  {
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const libMesh::Real rho = cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp];
    const libMesh::Real R_mix = cache.get_cached_values(Cache::MIXTURE_GAS_CONSTANT)[qp];
    const std::vector<libMesh::Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];

    this->check_and_reset_temp_cache(T);

    _kinetics.omega_dot( *(_temp_cache.get()), rho, R_mix, Y, omega_dot );

    return;
  }

  template<typename AntiochThermo>
  void AntiochEvaluator<AntiochThermo>::check_and_reset_temp_cache( const libMesh::Real T )
  {
    if( _temp_cache->T != T )
      {
        _temp_cache.reset( new Antioch::TempCache<libMesh::Real>(T) );
      }

    return;
  }

} // end namespace GRINS

#endif //GRINS_HAVE_ANTIOCH
