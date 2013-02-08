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

#ifndef GRINS_KINETICS_H
#define GRINS_KINETICS_H

// libMesh
#include "libmesh/getpot.h"

// GRINS
#include "grins/cached_values.h"
#include "grins/chemical_mixture.h"
#include "grins/reaction_set.h"

namespace GRINS
{

  class Kinetics
  {
  public:

    Kinetics( const GetPot& input, const ChemicalMixture& chem_mixture );
    ~Kinetics();
    
    //! Compute source terms
    void omega_dot( const CachedValues& cache, unsigned int qp,
		    std::vector<libMesh::Real>& omega_dot ) const;

    const ReactionSet& reaction_set() const;

  protected:

    void read_reaction_set_data_xml( const std::string& chem_file,
				     const bool verbose,
				     const ChemicalMixture& chem_mixture,
				     ReactionSet& reaction_set );

    ReactionSet _reaction_set;

  };

  /* ------------------------- Inline Functions -------------------------*/

  inline
  void Kinetics::omega_dot( const CachedValues& cache, unsigned int qp,
			    std::vector<libMesh::Real>& omega_dot ) const
  {
    _reaction_set.compute_mass_sources( cache.get_cached_values(Cache::TEMPERATURE)[qp],
					cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp],
					cache.get_cached_values(Cache::MIXTURE_GAS_CONSTANT)[qp],
					cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp],
					cache.get_cached_vector_values(Cache::MOLAR_DENSITIES)[qp],
					cache.get_cached_vector_values(Cache::SPECIES_NORMALIZED_ENTHALPY_MINUS_NORMALIZED_ENTROPY)[qp],
					omega_dot );
    return;
  }

  inline
  const ReactionSet& Kinetics::reaction_set() const
  {
    return _reaction_set;
  }

} // namespace GRINS

#endif // GRINS_KINETICS_H
