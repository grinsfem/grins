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

#include "grins_config.h"

#ifdef GRINS_HAVE_CANTERA

// C++
#include <vector>

// libMesh
#include "libmesh/libmesh_common.h"

// libMesh forward declarations
class GetPot;

// Cantera forward declarations
namespace Cantera
{
  class IdealGasMix;
  class Transport;
}

namespace GRINS
{

  // GRINS forward declarations
  class CachedValues;
  class CanteraMixture;

  class CanteraTransport
  {
  public:
    
    CanteraTransport( CanteraMixture& mixture );
    ~CanteraTransport();

    libMesh::Real mu( const CachedValues& cache, unsigned int qp ) const;

    libMesh::Real k( const CachedValues& cache, unsigned int qp ) const;

    void D( const CachedValues& cache, unsigned int qp,
	    std::vector<libMesh::Real>& D ) const;

  protected:

    Cantera::IdealGasMix& _cantera_gas;

    Cantera::Transport& _cantera_transport;

  private:

    CanteraTransport();

  };

} // namespace GRINS

#endif // GRINS_HAVE_CANTERA

#endif // GRINS_CANTERA_TRANSPORT_H
