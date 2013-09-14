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


#ifndef GRINS_ANTIOCH_EVALUATOR_H
#define GRINS_ANTIOCH_EVALUATOR_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// GRINS
#include "grins/antioch_mixture.h"
#include "grins/antioch_kinetics.h"
#include "grins/cached_values.h"
#include "grins/property_types.h"

// Antioch
#include "antioch/temp_cache.h"
#include "antioch/cea_evaluator.h"
#include "antioch/stat_mech_thermo.h"

// Boost
#include <boost/scoped_ptr.hpp>

namespace GRINS
{
  template<typename Thermo>
  class AntiochEvaluator
  {
  public:

    AntiochEvaluator( const AntiochMixture& mixture );

    virtual ~AntiochEvaluator();

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

    libMesh::Real cp( const libMesh::Real& T,
                      const std::vector<libMesh::Real>& Y );

    // Kinetics
    void omega_dot( const CachedValues& cache, unsigned int qp,
		    std::vector<libMesh::Real>& omega_dot );

  protected:

    const AntiochMixture& _chem;
    
    // This is a template type
    boost::scoped_ptr<Thermo> _thermo;

    boost::scoped_ptr<AntiochKinetics> _kinetics;

    boost::scoped_ptr<Antioch::TempCache<libMesh::Real> > _temp_cache;

    //! Helper method for managing _temp_cache
    /*! T *MUST* be pass-by-reference because of the structure
        of Antioch::TempCache! */
    void check_and_reset_temp_cache( const libMesh::Real& T );

    /* Below we will specialize the specialized_build_* functions to the appropriate type.
       This way, we can control how the cached transport objects get constructed
       based on the template type. This is achieved by the dummy types forcing operator
       overloading for each of the specialized types. */
    void build_thermo( const AntiochMixture& mixture )
    { specialized_build_thermo( mixture, _thermo, thermo_type<Thermo>() ); }

  private:

    AntiochEvaluator();

    void specialized_build_thermo( const AntiochMixture& mixture,
                                   boost::scoped_ptr<Antioch::StatMechThermodynamics<libMesh::Real> >& thermo,
                                   thermo_type<Antioch::StatMechThermodynamics<libMesh::Real> > )
    {
      thermo.reset( new Antioch::StatMechThermodynamics<libMesh::Real>( mixture.chemical_mixture() ) );
      return;
    }
    
    void specialized_build_thermo( const AntiochMixture& mixture,
                                   boost::scoped_ptr<Antioch::CEAEvaluator<libMesh::Real> >& thermo,
                                   thermo_type<Antioch::CEAEvaluator<libMesh::Real> > )
    {
      thermo.reset( new Antioch::CEAEvaluator<libMesh::Real>( mixture.cea_mixture() ) );
      return;
    }

  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename Thermo>
  inline
  libMesh::Real AntiochEvaluator<Thermo>::M( unsigned int species ) const
  {
    return _chem.M(species);
  }

  template<typename Thermo>
  inline
  libMesh::Real AntiochEvaluator<Thermo>::M_mix( const std::vector<libMesh::Real>& mass_fractions ) const
  {
    return _chem.M_mix(mass_fractions);
  }

  template<typename Thermo>
  inline
  libMesh::Real AntiochEvaluator<Thermo>::R( unsigned int species ) const
  {
    return _chem.R(species);
  }

  template<typename Thermo>
  inline
  libMesh::Real AntiochEvaluator<Thermo>::R_mix( const std::vector<libMesh::Real>& mass_fractions ) const
  {
    return _chem.R_mix(mass_fractions);
  }
  
  template<typename Thermo>
  inline
  libMesh::Real AntiochEvaluator<Thermo>::X( unsigned int species, libMesh::Real M, libMesh::Real mass_fraction ) const
  {
    return _chem.X(species,M,mass_fraction);
  }
  
  template<typename Thermo>
  inline
  void AntiochEvaluator<Thermo>::X( libMesh::Real M, const std::vector<libMesh::Real>& mass_fractions, 
                                    std::vector<libMesh::Real>& mole_fractions ) const
  {
    _chem.X(M,mass_fractions,mole_fractions);
    return;
  }
  
  template<typename Thermo>
  inline
  unsigned int AntiochEvaluator<Thermo>::species_index( const std::string& species_name ) const
  {
    return _chem.species_index(species_name);
  }
  
  template<typename Thermo>
  inline
  std::string AntiochEvaluator<Thermo>::species_name( unsigned int species_index ) const
  {
    return _chem.species_name(species_index);
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_EVALUATOR_H
