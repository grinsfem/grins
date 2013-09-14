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

// This class
#include "grins/antioch_wilke_transport_evaluator.h"

// GRINS
#include "grins/cached_values.h"

namespace GRINS
{
  template<typename Thermo, typename Viscosity, typename Conductivity, typename Diffusivity>
  AntiochWilkeTransportEvaluator<Thermo,Viscosity,Conductivity,Diffusivity>::AntiochWilkeTransportEvaluator( const AntiochWilkeTransportMixture<Thermo,Viscosity,Conductivity,Diffusivity>& mixture )
    : AntiochEvaluator<Thermo>( mixture ),
      _wilke_evaluator( new Antioch::WilkeEvaluator<Viscosity,Conductivity>( mixture.wilke_mixture(), mixture.viscosity(), mixture.conductivity() ) ),
      _diffusivity( mixture.diffusivity() )
  {
    return;
  }

  template<typename T, typename V, typename C, typename D>
  AntiochWilkeTransportEvaluator<T,V,C,D>::~AntiochWilkeTransportEvaluator()
  {
    return;
  }

  template<typename Th, typename V, typename C, typename D>
  libMesh::Real AntiochWilkeTransportEvaluator<Th,V,C,D>::mu( const CachedValues& cache, unsigned int qp )
  {
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const std::vector<libMesh::Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];

    return this->mu( T, Y );
  }

  template<typename Th, typename V, typename C, typename D>
  libMesh::Real AntiochWilkeTransportEvaluator<Th,V,C,D>::k( const CachedValues& cache, unsigned int qp )
  {
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const std::vector<libMesh::Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];

    return this->k( T, Y );
  }

  template<typename Th, typename V, typename C, typename D>
  void AntiochWilkeTransportEvaluator<Th,V,C,D>::mu_and_k( const CachedValues& cache, unsigned int qp,
                                                           libMesh::Real& mu, libMesh::Real& k ) 
  {
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const std::vector<libMesh::Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];

    _wilke_evaluator->mu_and_k( T, Y, mu, k );
    return;
  }

  template<typename Th, typename V, typename C, typename D>
  void AntiochWilkeTransportEvaluator<Th,V,C,D>::D( const CachedValues& cache, unsigned int qp,
                                                    std::vector<libMesh::Real>& D )
  {
    const libMesh::Real rho = cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp];
    
    /*! \todo Find a way to cache these so we don't have to recompute them */
    const libMesh::Real cp = this->cp(cache,qp);
    const libMesh::Real k = this->k(cache,qp);

    this->D(rho,cp,k,D);
    
    return;
  }

  template<typename Th, typename V, typename C, typename D>
  libMesh::Real AntiochWilkeTransportEvaluator<Th,V,C,D>::mu( const libMesh::Real T,
                                                              const std::vector<libMesh::Real>& Y )
  {
    return _wilke_evaluator->mu( T, Y );
  }
  
  template<typename Th, typename V, typename C, typename D>
  libMesh::Real AntiochWilkeTransportEvaluator<Th,V,C,D>::k( const libMesh::Real T,
                                                             const std::vector<libMesh::Real>& Y )
  {
    return _wilke_evaluator->k( T, Y );
  }

  template<typename Th, typename V, typename C, typename D>
  void AntiochWilkeTransportEvaluator<Th,V,C,D>::D( const libMesh::Real rho, const libMesh::Real cp,
                                                    const libMesh::Real k,
                                                    std::vector<libMesh::Real>& D )
  {
    std::fill( D.begin(), D.end(), _diffusivity.D(rho,cp,k) );

    return;
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
