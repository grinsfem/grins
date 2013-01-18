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

#ifndef GRINS_CEA_THERMO_H
#define GRINS_CEA_THERMO_H

// C++
#include <iomanip>

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/libmesh_common.h"

// GRINS
#include "cached_values.h"
#include "chemical_mixture.h"
#include "cea_curve_fit.h"

namespace GRINS
{
  class CEAThermodynamics
  {
  public:

    CEAThermodynamics( const GetPot& input, const ChemicalMixture& chem_mixture );
    ~CEAThermodynamics();

    Real cp( Real T, unsigned int species ) const;

    Real cp( Real T, const std::vector<Real>& mass_fractions ) const;

    Real cv( Real T, unsigned int species ) const;

    Real cv( Real T, const std::vector<Real>& mass_fractions ) const;

    Real h( Real T, unsigned int species ) const;

    void h( Real T, std::vector<Real>& h ) const;

    Real h_RT_minus_s_R( Real T, unsigned int species ) const;

    void h_RT_minus_s_R( Real T, std::vector<Real>& h_RT_minus_s_R ) const;

    Real cp_over_R( Real T, unsigned int species ) const;

    Real h_over_RT( Real T, unsigned int species ) const;

    Real s_over_R( Real T, unsigned int species ) const;

    /* -------------- Ideal Gas Mixture Interaface Methods --------------*/

    Real cp( const CachedValues& cache, unsigned int qp ) const;

    Real cv( const CachedValues& cache, unsigned int qp ) const;
     
    Real h(const CachedValues& cache, unsigned int qp, unsigned int species) const;

    void h(const CachedValues& cache, unsigned int qp, std::vector<Real>& h) const;

    Real h_RT_minus_s_R( const CachedValues& cache, unsigned int qp,
			 unsigned int species ) const;

    void h_RT_minus_s_R( const CachedValues& cache, unsigned int qp,
			 std::vector<Real>& h_RT_minus_s_R) const;

  protected:

    void read_thermodynamic_table();

    void read_thermodynamic_table( std::istream& in );

    const ChemicalMixture& _chem_mixture;

    std::vector<CEACurveFit*> _species_curve_fits;

    std::vector<Real> _cp_at_200p1;

  private:
    
    CEAThermodynamics();

  };
  
  /* ------------------------- Inline Functions -------------------------*/
  
  inline
  Real CEAThermodynamics::cv( Real T, unsigned int species ) const
  { return this->cp(T,species) - _chem_mixture.R(species); }

  inline
  Real CEAThermodynamics::cv( Real T, const std::vector<Real>& mass_fractions ) const
  { return this->cp(T,mass_fractions) - _chem_mixture.R(mass_fractions); }

  inline
  Real CEAThermodynamics::cp( const CachedValues& cache, unsigned int qp ) const
  {
    return this->cp( cache.get_cached_values(Cache::TEMPERATURE)[qp],
		     cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp] );
  }

  inline
  Real CEAThermodynamics::cv( const CachedValues& cache, unsigned int qp ) const
  {
    return this->cv( cache.get_cached_values(Cache::TEMPERATURE)[qp],
		     cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp] );
  }
     
  inline
  Real CEAThermodynamics::h(const CachedValues& cache, unsigned int qp, 
			    unsigned int species) const
  {
    return this->h( cache.get_cached_values(Cache::TEMPERATURE)[qp], species );
  }

  inline
  void CEAThermodynamics::h(const CachedValues& cache, unsigned int qp,
			    std::vector<Real>& h) const
  {
    this->h( cache.get_cached_values(Cache::TEMPERATURE)[qp], h );
    return;
  }

  inline
  Real CEAThermodynamics::h_RT_minus_s_R( const CachedValues& cache, unsigned int qp,
					  unsigned int species ) const
  {
    return this->h_RT_minus_s_R( cache.get_cached_values(Cache::TEMPERATURE)[qp],
				 species );
  }

  inline
  void CEAThermodynamics::h_RT_minus_s_R( const CachedValues& cache, unsigned int qp,
					  std::vector<Real>& h_RT_minus_s_R) const
  {
    this->h_RT_minus_s_R( cache.get_cached_values(Cache::TEMPERATURE)[qp],
			  h_RT_minus_s_R );
    return;
  }
  

} // namespace GRINS

#endif //GRINS_CEA_THERMO_H
