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

#ifndef GRINS_ANTIOCH_EVALUATOR_H
#define GRINS_ANTIOCH_EVALUATOR_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// GRINS
#include "grins/antioch_mixture.h"
#include "grins/antioch_kinetics.h"

// Antioch
#include "antioch/temp_cache.h"

// Boost
#include <boost/scoped_ptr.hpp>

namespace GRINS
{
  // GRINS forward declarations
  class CachedValues;

  template<typename Thermo, typename Transport>
  class AntiochEvaluator
  {
  public:

    AntiochEvaluator( const AntiochMixture& mixture );
    ~AntiochEvaluator();

    // Chemistry
    libMesh::Real M( unsigned int species ) const;

    libMesh::Real M_mix( const std::vector<libMesh::Real>& mass_fractions ) const;

    libMesh::Real R( unsigned int species ) const;

    libMesh::Real R_mix( const std::vector<libMesh::Real>& mass_fractions ) const;

    libMesh::Real X( unsigned int species, libMesh::Real M, libMesh::Real mass_fraction ) const;

    void X( libMesh::Real M, const std::vector<libMesh::Real>& mass_fractions, 
	    std::vector<libMesh::Real>& mole_fractions ) const;

    unsigned int species_index( const std::string& species_name ) const;

    std::string species_name( unsigned int species_index ) const;

    // Thermo
    libMesh::Real cp( const CachedValues& cache, unsigned int qp );

    libMesh::Real cv( const CachedValues& cache, unsigned int qp );
     
    libMesh::Real h_s(const CachedValues& cache, unsigned int qp, unsigned int species);

    void h_s(const CachedValues& cache, unsigned int qp, std::vector<libMesh::Real>& h_s);

    // Transport
    libMesh::Real mu( const CachedValues& cache, unsigned int qp );

    libMesh::Real k( const CachedValues& cache, unsigned int qp );

    void mu_and_k( const CachedValues& cache, unsigned int qp,
                   libMesh::Real& mu, libMesh::Real k );

    void D( const CachedValues& cache, unsigned int qp,
	    std::vector<libMesh::Real>& D );

    // Kinetics
    void omega_dot( const CachedValues& cache, unsigned int qp,
		    std::vector<libMesh::Real>& omega_dot );

  protected:

    const AntiochMixture& _chem;
    
    // This is a template type
    Thermo _thermo;

    // This is a template type
    Transport _transport;

    AntiochKinetics _kinetics;

    boost::scoped_ptr<Antioch::TempCache<libMesh::Real> > _temp_cache;

    //! Helper method for managing _temp_cache
    /*! T *MUST* be pass-by-reference because of the structure
        of Antioch::TempCache! */
    void check_and_reset_temp_cache( const libMesh::Real& T );

  private:

    AntiochEvaluator();

  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename Thermo, typename Transport>
  inline
  libMesh::Real AntiochEvaluator<Thermo,Transport>::M( unsigned int species ) const
  {
    return _chem.M(species);
  }

  template<typename Thermo, typename Transport>
  inline
  libMesh::Real AntiochEvaluator<Thermo,Transport>::M_mix( const std::vector<libMesh::Real>& mass_fractions ) const
  {
    return _chem.M_mix(mass_fractions);
  }

  template<typename Thermo, typename Transport>
  inline
  libMesh::Real AntiochEvaluator<Thermo,Transport>::R( unsigned int species ) const
  {
    return _chem.R(species);
  }

  template<typename Thermo, typename Transport>
  inline
  libMesh::Real AntiochEvaluator<Thermo,Transport>::R_mix( const std::vector<libMesh::Real>& mass_fractions ) const
  {
    return _chem.R_mix(mass_fractions);
  }
  
  template<typename Thermo, typename Transport>
  inline
  libMesh::Real AntiochEvaluator<Thermo,Transport>::X( unsigned int species, libMesh::Real M, libMesh::Real mass_fraction ) const
  {
    return _chem.X(species,M,mass_fraction);
  }
  
  template<typename Thermo, typename Transport>
  inline
  void AntiochEvaluator<Thermo,Transport>::X( libMesh::Real M, const std::vector<libMesh::Real>& mass_fractions, 
                            std::vector<libMesh::Real>& mole_fractions ) const
  {
    _chem.X(M,mass_fractions,mole_fractions);
    return;
  }
  
  template<typename Thermo, typename Transport>
  inline
  unsigned int AntiochEvaluator<Thermo,Transport>::species_index( const std::string& species_name ) const
  {
    return _chem.species_index(species_name);
  }
  
  template<typename Thermo, typename Transport>
  inline
  std::string AntiochEvaluator<Thermo,Transport>::species_name( unsigned int species_index ) const
  {
    return _chem.species_name(species_index);
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_EVALUATOR_H
