//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_CONSTANT_TRANSPORT_H
#define GRINS_CONSTANT_TRANSPORT_H

// libMesh
#include "getpot.h"

// GRINS
#include "chemical_mixture.h"
#include "reacting_flow_cache.h"
#include "cached_values.h"

namespace GRINS
{
  class ConstantTransport
  {
  public:
    
    ConstantTransport( const GetPot& input, const ChemicalMixture& chem_mixture );
    ~ConstantTransport();

    inline
    Real mu( const ReactingFlowCache& )
    { return _mu; }

    inline
    Real k( const ReactingFlowCache& )
    { return _k; }
    
    inline
    void D( const ReactingFlowCache& cache, std::vector<Real>& D )
    { std::fill( D.begin(), D.end(), _Le*_k/( cache.rho() * cache.cp() ) ); }

    Real mu( const CachedValues& cache, unsigned int qp );

    Real k( const CachedValues& cache, unsigned int qp );

    void D( const CachedValues& cache, unsigned int qp,
	    std::vector<Real>& D );

  protected:
    
    const ChemicalMixture& _chem_mixture;

    const Real _mu;
    const Real _k;
    const Real _Le; 

  private:
    
    ConstantTransport();

  };

  inline
  Real ConstantTransport::mu( const CachedValues& /*cache*/, unsigned int /*qp*/ )
  {
    return _mu;
  }

  inline
  Real ConstantTransport::k( const CachedValues& /*cache*/, unsigned int /*qp*/ )
  {
    return _k;
  }

  inline
  void ConstantTransport::D( const CachedValues& cache, unsigned int qp,
			     std::vector<Real>& D )
  { const Real rho = cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp];
    const Real cp  = cache.get_cached_values(Cache::MIXTURE_SPECIFIC_HEAT_P)[qp];
    std::fill( D.begin(), D.end(), _Le*_k/( rho*cp ) );
    return;
  }
  
} // namespace GRINS

#endif //GRINS_CONSTANT_TRANSPORT_H
