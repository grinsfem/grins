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

#ifndef GRINS_ANTIOCH_STAT_MECH_THERMO_H
#define GRINS_ANTIOCH_STAT_MECH_THERMO_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// GRINS
#include "grins/antioch_mixture.h"

// libMesh
#include "libmesh/libmesh_common.h"

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/vector_utils.h"
#include "antioch/stat_mech_thermo.h"
#include "antioch/temp_cache.h"

namespace GRINS
{
  class AntiochStatMechThermo
  {
  public:

    AntiochStatMechThermo( const AntiochMixture& mixture );
    ~AntiochStatMechThermo();

    libMesh::Real cp( const Antioch::TempCache<libMesh::Real>& cache, unsigned int species ) const;

    libMesh::Real cp( const Antioch::TempCache<libMesh::Real>& cache,
                      const std::vector<libMesh::Real> mass_fractions ) const;

    libMesh::Real cv( const Antioch::TempCache<libMesh::Real>& cache, unsigned int species ) const;

    libMesh::Real cv( const Antioch::TempCache<libMesh::Real>& cache,
                      const std::vector<libMesh::Real> mass_fractions ) const;

    libMesh::Real h_s( const Antioch::TempCache<libMesh::Real>& cache, unsigned int species ) const;

    void h_s( const Antioch::TempCache<libMesh::Real>& cache, std::vector<libMesh::Real>& h_s ) const;

  protected:

    const AntiochMixture& _antioch_mixture;

    Antioch::StatMechThermodynamics<libMesh::Real> _antioch_thermo;

  private:

    AntiochStatMechThermo();
    
  };

  /* ------------------------- Inline Functions -------------------------*/
  inline
  libMesh::Real AntiochStatMechThermo::cp( const Antioch::TempCache<libMesh::Real>& /*cache*/,
                                           unsigned int /*species*/ ) const
  {
    // Not yet implemented in Antioch
    libmesh_not_implemented();
    return 0.0;
  }

  inline
  libMesh::Real AntiochStatMechThermo::cp( const Antioch::TempCache<libMesh::Real>& cache,
                                           const std::vector<libMesh::Real> mass_fractions ) const
  {
    return _antioch_thermo.cp( cache.T, cache.T, mass_fractions );
  }

  inline
  libMesh::Real AntiochStatMechThermo::cv( const Antioch::TempCache<libMesh::Real>& cache, unsigned int species ) const
  {
    return _antioch_thermo.cv( species, cache.T, cache.T );
  }

  inline
  libMesh::Real AntiochStatMechThermo::cv( const Antioch::TempCache<libMesh::Real>& cache,
                                      const std::vector<libMesh::Real> mass_fractions ) const
  {
    return _antioch_thermo.cv( cache.T, cache.T, mass_fractions );
  }

  inline
  libMesh::Real AntiochStatMechThermo::h_s( const Antioch::TempCache<libMesh::Real>& cache,
                                            unsigned int species ) const
  {
    return _antioch_thermo.h_tot( species, cache.T);
  }

  inline
  void AntiochStatMechThermo::h_s( const Antioch::TempCache<libMesh::Real>& cache,
                                   std::vector<libMesh::Real>& h_s ) const
  {
    for( unsigned int s = 0; s < _antioch_mixture.n_species(); s++ )
      {
        h_s[s] = this->h_s(cache,s);
      }

    return;
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_STAT_MECH_THERMO_H
