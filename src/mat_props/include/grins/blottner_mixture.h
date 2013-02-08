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

#ifndef GRINS_BLOTTNER_MIXTURE_H
#define GRINS_BLOTTNER_MIXTURE_H

// GRINS
#include "grins/blottner_viscosity.h"

namespace GRINS
{

  class BlottnerMixture
  {
  public:

    BlottnerMixture( const GetPot& input, const ChemicalMixture& chem_mixture );
    ~BlottnerMixture();

    libMesh::Real mu( const CachedValues& cache, unsigned int qp,
	     unsigned int species ) const;

    libMesh::Real mu( libMesh::Real T, unsigned int species ) const;
    libMesh::Real mu( libMesh::Real T ) const;

  protected:

    void read_blottner_table();

    const ChemicalMixture& _chem_mixture;

    std::vector<BlottnerViscosity*> _species_viscosities;

  };

  inline
  libMesh::Real BlottnerMixture::mu( const CachedValues& cache, unsigned int qp,
			    unsigned int species ) const
  { return this->mu( cache.get_cached_values(Cache::TEMPERATURE)[qp], species ) };

} // namespace GRINS

#endif // GRINS_BLOTTNER_MIXTURE_H
