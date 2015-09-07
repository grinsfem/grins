//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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
#include "grins/antioch_constant_transport_evaluator.h"

// GRINS
#include "grins/cached_values.h"

namespace GRINS
{
  template<typename Thermo, typename Conductivity>
  AntiochConstantTransportEvaluator<Thermo,Conductivity>::AntiochConstantTransportEvaluator( const AntiochConstantTransportMixture<Conductivity>& mixture )
    : AntiochEvaluator<Thermo>( mixture ),
      _mu( mixture.mu() ),
      _conductivity( mixture.conductivity() ),
      _diffusivity( mixture.diffusivity() )
  {
    return;
  }

  template<typename Thermo, typename Conductivity>
  AntiochConstantTransportEvaluator<Thermo,Conductivity>::~AntiochConstantTransportEvaluator()
  {
    return;
  }

  template<typename Thermo, typename Conductivity>
  libMesh::Real AntiochConstantTransportEvaluator<Thermo,Conductivity>::mu( const CachedValues& /*cache*/,
                                                                            unsigned int /*qp*/ )
  {
    return this->_mu;
  }

  template<typename Thermo, typename Conductivity>
  libMesh::Real AntiochConstantTransportEvaluator<Thermo,Conductivity>::k( const CachedValues& cache,
                                                                           unsigned int qp )
  {
    libMesh::Real cp = this->cp( cache, qp );

    return _conductivity( _mu, cp );
  }

  template<typename Thermo, typename Conductivity>
  void AntiochConstantTransportEvaluator<Thermo,Conductivity>::mu_and_k( const CachedValues& cache,
                                                                         unsigned int qp,
                                                                         libMesh::Real& mu,
                                                                         libMesh::Real& k ) 
  {
    mu = _mu;

    libMesh::Real cp = this->cp( cache, qp );

    k = _conductivity( _mu, cp );
    return;
  }

  template<typename Thermo, typename Conductivity>
  void AntiochConstantTransportEvaluator<Thermo,Conductivity>::D( const CachedValues& cache,
                                                                  unsigned int qp,
                                                                  std::vector<libMesh::Real>& D )
  {
    const libMesh::Real rho = cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp];
    
    /*! \todo Find a way to cache these so we don't have to recompute them */
    const libMesh::Real cp = this->cp(cache,qp);
    
    const libMesh::Real k = _conductivity( _mu, cp );

    this->D(rho,cp,k,D);
    
    return;
  }

  template<typename Thermo, typename Conductivity>
  libMesh::Real AntiochConstantTransportEvaluator<Thermo,Conductivity>::mu( const libMesh::Real /*T*/,
                                                                            const std::vector<libMesh::Real>& /*Y*/ )
  {
    return _mu;
  }
  
  template<typename Thermo, typename Conductivity>
  libMesh::Real AntiochConstantTransportEvaluator<Thermo,Conductivity>::k( const libMesh::Real T,
                                                                           const std::vector<libMesh::Real>& Y )
  {
    const libMesh::Real cp = this->cp( T, Y );

    return _conductivity( _mu, cp );
  }

  template<typename Thermo, typename Conductivity>
  void AntiochConstantTransportEvaluator<Thermo,Conductivity>::D( const libMesh::Real rho, const libMesh::Real cp,
                                                                  const libMesh::Real k,
                                                                  std::vector<libMesh::Real>& D )
  {
    std::fill( D.begin(), D.end(), _diffusivity.D(rho,cp,k) );
  }

  template<typename Thermo, typename Conductivity>
  void AntiochConstantTransportEvaluator<Thermo,Conductivity>::mu_and_k_and_D( const libMesh::Real /*T*/,
                                                                               const libMesh::Real rho,
                                                                               const libMesh::Real cp,
                                                                               const std::vector<libMesh::Real>& /*Y*/,
                                                                               libMesh::Real& mu, libMesh::Real& k,
                                                                               std::vector<libMesh::Real>& D )
  {
    mu = _mu;

    k = _conductivity( _mu, cp );

    std::fill( D.begin(), D.end(), _diffusivity.D(rho,cp,k) );
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
