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

#ifndef GRINS_CONSTANT_TRANSPORT_H
#define GRINS_CONSTANT_TRANSPORT_H

// libMesh
#include "libmesh/getpot.h"

// GRINS
#include "grins/chemical_mixture.h"
#include "grins/cached_values.h"

namespace GRINS
{
  class ConstantTransport
  {
  public:
    
    ConstantTransport( const GetPot& input, const ChemicalMixture& chem_mixture );
    ~ConstantTransport();

    libMesh::Real mu( const CachedValues& cache, unsigned int qp ) const;

    libMesh::Real k( const CachedValues& cache, unsigned int qp ) const;

    void D( const CachedValues& cache, unsigned int qp,
	    std::vector<libMesh::Real>& D ) const;

  protected:
    
    const ChemicalMixture& _chem_mixture;

    const libMesh::Real _mu;
    const libMesh::Real _k;
    const libMesh::Real _Le; 

  private:
    
    ConstantTransport();

  };

  inline
  libMesh::Real ConstantTransport::mu( const CachedValues& /*cache*/, unsigned int /*qp*/ ) const
  {
    return _mu;
  }

  inline
  libMesh::Real ConstantTransport::k( const CachedValues& /*cache*/, unsigned int /*qp*/ ) const
  {
    return _k;
  }

  inline
  void ConstantTransport::D( const CachedValues& cache, unsigned int qp,
			     std::vector<libMesh::Real>& D ) const
  { const libMesh::Real rho = cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp];
    const libMesh::Real cp  = cache.get_cached_values(Cache::MIXTURE_SPECIFIC_HEAT_P)[qp];
    std::fill( D.begin(), D.end(), _Le*_k/( rho*cp ) );
    return;
  }
  
} // namespace GRINS

#endif //GRINS_CONSTANT_TRANSPORT_H
