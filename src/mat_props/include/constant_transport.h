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

  protected:
    
    const ChemicalMixture& _chem_mixture;

    const Real _mu;
    const Real _k;
    const Real _Le; 

  private:
    
    ConstantTransport();

  };

} // namespace GRINS

#endif //GRINS_CONSTANT_TRANSPORT_H
