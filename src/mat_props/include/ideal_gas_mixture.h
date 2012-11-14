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

// GRINS
#include "reacting_flow_cache.h"
#include "chemical_mixture.h"

// Thermo classes
#include "cea_thermo.h"
#include "cantera_thermo.h"

// Transport classes
#include "cantera_transport.h"

// Kinetics classes
#include "cantera_kinetics.h"

namespace GRINS
{
  template<typename Thermo, typename Transport, typename Kinetics>
  class IdealGasMixture
  {
  public:

    IdealGasMixture( const GetPot& input );
    ~IdealGasMixture();

    // Chemistry quantities
    inline
    Real R( unsigned int species ) const
    { return this->_chem_mixture.R(species); }

    inline
    Real R( const std::vector<Real>& mass_fractions ) const
    { return this->_chem_mixture.R(mass_fractions); }

    inline
    Real M( unsigned int species ) const
    { return this->_chem_mixture.M(species); }

    inline
    Real M( const std::vector<Real>& mass_fractions ) const
    { return this->_chem_mixture.M(mass_fractions); }

    inline
    Real X( unsigned int species, Real M, Real mass_fraction ) const
    { return this->_chem_mixture.X(species,M,mass_fraction); }

    inline
    void X( Real M, const std::vector<Real>& mass_fractions, 
	    std::vector<Real>& mole_fractions ) const
    { return this->_chem_mixture.X(M,mass_fractions,mole_fractions); }


    // Thermodynamic quantities
    inline
    Real cp( const ReactingFlowCache& cache )
    { return this->_thermo.cp(cache); }

    inline
    Real cv( const ReactingFlowCache& cache )
    { return this->_thermo.cv(cache); }

    // Transport quantities
    inline
    Real mu( const ReactingFlowCache& cache )
    { return this->_transport.mu(cache); }

    inline
    Real k( const ReactingFlowCache& cache )
    { return this->_transport.k(cache); }

    inline
    void D( const ReactingFlowCache& cache, std::vector<Real>& D )
    { return this->_transport.D(cache,D); }

    // Kinetics quantites
    void omega_dot( const ReactingFlowCache& cache, std::vector<Real>& omega_dot )
    { return this->_kinetics.omega_dot(cache,omega_dot); }

  protected:
    
    std::vector<std::string> read_species_list( const GetPot& input );

    //! Stash data for the requested chemical mixture
    /*!
     * A reference will be cached in Thermo, Transport, and Kinetics objects as
     * needed. The intent is for this object to maintain ownership.
     */
    ChemicalMixture _chem_mixture;
    
    Thermo _thermo;
    
    Transport _transport;

    Kinetics _kinetics;

  };
}

#endif //GRINS_IDEAL_GAS_MIXTURE_H
