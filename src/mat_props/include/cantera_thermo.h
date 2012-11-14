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

#ifndef GRINS_CANTERA_THERMO_H
#define GRINS_CANTERA_THERMO_H

// libMesh
#include "getpot.h"

// GRINS
#include "chemical_mixture.h"
#include "reacting_flow_cache.h"
#include "cantera_singleton.h"

#ifdef HAVE_CANTERA

namespace GRINS
{

  class CanteraThermodynamics
  {
  public:

    CanteraThermodynamics( const GetPot& input, const ChemicalMixture& chem_mixture );
    ~CanteraThermodynamics();

    Real cp( const ReactingFlowCache& cache );

    Real cv( const ReactingFlowCache& cache );

    Real h(const ReactingFlowCache& cache, unsigned int species);

    void h(const ReactingFlowCache& cache, std::vector<Real>& h);

  protected:

    const ChemicalMixture& _chem_mixture;

    Cantera::IdealGasMix& _cantera_gas;

  private:

    CanteraThermodynamics();

  };

}

#endif //HAVE_CANTERA

#endif //GRINS_CANTERA_THERMO_H
