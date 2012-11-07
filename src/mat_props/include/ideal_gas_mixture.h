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

#ifndef GRINS_IDEAL_GAS_MIXTURE_H
#define GRINS_IDEAL_GAS_MIXTURE_H

namespace GRINS
{
  template<typename Thermo, typename Transport, typename Kinetics>
  class IdealGasMixture
  {
  public:

    IdealGasMixture( const GetPot& input );
    ~IdealGasMixture();

    // Thermodynamic quantities
    void cp( const CachedQuantities& cache, std::vector<Real>& cp_values );
    Real cp( const CachedQuantities& cache, unsigned int species );

    void cv( const CachedQuantities& cache, std::vector<Real>& cp_values );
    Real cv( const CachedQuantities& cache, unsigned int species );

    Real R( unsigned int species );
    Real R( const std::vector<Real>& mass_fractions );

    Real M( unsigned int species );
    Real M( const std::vector<Real>& mass_fractions );

    Real X( Real M, Real mass_fraction  );
    void X( Real M, const std::vector<Real>& mass_fractions, std::vector<Real>& mole_fractions );

    // Transport quantities
    void mu( const CachedQuantities& cache, std::vector<Real>& mu_values );
    Real mu( const CachedQuantities& cache, unsigned int species );

    void k( const CachedQuantities& cache, std::vector<Real>& k_values );
    Real k( const CachedQuantities& cache, unsigned int species );

    // Kinetics quantites
    void omega_dot( const CachedQuantities& cache, std::vector<Real>& omega );
    Real omega_dot( const CachedQuantities& cache, unsigned int species );

  protected:
    
    std::vector<std::string> read_species_list( const GetPot& input );

    Thermo _thermo;
    
    Transport _transport;

    Kinetics _kinetics;
    
    //! Stash data for the requested chemical mixture
    /*!
     * A pointer will be passed to Thermo, Transport, and Kinetics objects as
     * needed.
     */
    std::tr1::shared_ptr<ChemicalMixture> _chem_mixture;

  };
}

#endif //GRINS_IDEAL_GAS_MIXTURE_H
