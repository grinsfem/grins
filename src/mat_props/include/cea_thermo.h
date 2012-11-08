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

#ifndef GRINS_CEA_THERMO_H
#define GRINS_CEA_THERMO_H

// C++
#include <iomanip>

// libMesh
#include "getpot.h"
#include "libmesh_common.h"

// GRINS
#include "chemical_mixture.h"
#include "cea_curve_fit.h"
#include "reacting_flow_cache.h"

namespace GRINS
{
  class CEAThermodynamics
  {
  public:

    CEAThermodynamics( const GetPot& input, const ChemicalMixture& chem_mixture );
    ~CEAThermodynamics();

    inline
    Real cp( const ReactingFlowCache& cache, unsigned int species ) const
    { return this->cp(cache.T(),species); }

    inline
    Real cp( const ReactingFlowCache& cache ) const
    { return this->cp(cache.T(),cache.mass_fractions()); }

    Real cp( Real T, unsigned int species ) const;

    Real cp( Real T, const std::vector<Real>& mass_fractions ) const;

    inline
    Real cv( const ReactingFlowCache& cache, unsigned int species ) const
    { return this->cv(cache.T(),species); }

    inline
    Real cv( const ReactingFlowCache& cache ) const
    { return this->cv(cache.T(),cache.mass_fractions()); }

    inline
    Real cv( Real T, unsigned int species ) const
    { return this->cp(T,species) - _chem_mixture.R(species); }

    inline
    Real cv( Real T, const std::vector<Real>& mass_fractions ) const
    { return this->cp(T,mass_fractions) - _chem_mixture.R(mass_fractions); }

  protected:

    void read_thermodynamic_table();

    void read_thermodynamic_table( std::istream& in );

    Real cp_over_R( Real T, unsigned int species ) const;

    const ChemicalMixture& _chem_mixture;

    std::vector<CEACurveFit*> _species_curve_fits;

    std::vector<Real> _cp_at_200p1;

  private:
    
    CEAThermodynamics();

  };

} // namespace GRINS

#endif //GRINS_CEA_THERMO_H
