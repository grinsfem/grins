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

#ifndef GRINS_CANTERA_TRANSPORT_H
#define GRINS_CANTERA_TRANSPORT_H

// libMesh
#include "libmesh/getpot.h"

// GRINS
#include "grins/chemical_mixture.h"
#include "grins/cached_values.h"
#include "grins/cantera_singleton.h"

// Cantera
#ifdef GRINS_HAVE_CANTERA
// Needs to be included *after* cantera_singleton.h (i.e. cantera/IdealGasMix.h)
#include "cantera/transport.h"
#endif

#ifdef GRINS_HAVE_CANTERA

namespace GRINS
{

  class CanteraTransport
  {
  public:
    
    CanteraTransport( const GetPot& input, const ChemicalMixture& chem_mixture );
    ~CanteraTransport();

    libMesh::Real mu( const CachedValues& cache, unsigned int qp ) const;

    libMesh::Real k( const CachedValues& cache, unsigned int qp ) const;

    void D( const CachedValues& cache, unsigned int qp,
	    std::vector<libMesh::Real>& D ) const;

  protected:

    const ChemicalMixture& _chem_mixture;

    Cantera::IdealGasMix& _cantera_gas;

    Cantera::Transport* _cantera_transport;

  private:

    CanteraTransport();

  };

} // namespace GRINS

#endif // GRINS_HAVE_CANTERA

#endif // GRINS_CANTERA_TRANSPORT_H
