//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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
#include "antioch/ideal_gas_thermo.h"

namespace GRINS
{
  //! Wrapper class for evaluating chemistry and thermo properties using Antioch
  /*!
    This class is expected to be constructed *after* threads have been forked and will only
    live during the lifetime of the thread.
    By default, Antioch is working in SI units. Note that this documentation will always
    be built regardless if Antioch is included in the GRINS build or not. Check configure
    output to confirm that Antioch was included in the build.
  */
  template<typename KineticsThermoCurveFit, typename Thermo>
  class AntiochEvaluator
  {
  public:

    AntiochEvaluator() = delete;

    AntiochEvaluator( const AntiochMixture<KineticsThermoCurveFit>& mixture );

    virtual ~AntiochEvaluator() = default;

    // Chemistry
    libMesh::Real M( unsigned int species ) const;

    libMesh::Real M_mix( const std::vector<libMesh::Real>& mass_fractions ) const;

    libMesh::Real R( unsigned int species ) const;

    libMesh::Real R_mix( const std::vector<libMesh::Real>& mass_fractions ) const;

    libMesh::Real X( unsigned int species, libMesh::Real M, libMesh::Real mass_fraction ) const;

    void X( libMesh::Real M, const std::vector<libMesh::Real> & mass_fractions,
            std::vector<libMesh::Real> & mole_fractions ) const;

    unsigned int species_index( const std::string& species_name ) const;

    std::string species_name( unsigned int species_index ) const;

    // Thermo
    libMesh::Real cp( const libMesh::Real & T, const libMesh::Real P, const std::vector<libMesh::Real> & Y )
    { return this->specialized_cp(T,P,Y,*_thermo); }

    void cp_s( const libMesh::Real & T,
               const libMesh::Real P,
               const std::vector<libMesh::Real> & Y,
               std::vector<libMesh::Real> & cp_s )
    { this->specialized_cp_s(T,P,Y,*_thermo,cp_s); }

    libMesh::Real cv( const libMesh::Real& T, const libMesh::Real P, const std::vector<libMesh::Real>& Y )
    { return this->specialized_cv(T,P,Y,*_thermo); }

    libMesh::Real h_s( const libMesh::Real & T, unsigned int species )
    { return this->specialized_h_s(T,species,*_thermo); }

    // Kinetics
    void omega_dot( const libMesh::Real& T, libMesh::Real rho,
                    const std::vector<libMesh::Real> mass_fractions,
                    std::vector<libMesh::Real>& omega_dot );

  protected:

    const AntiochMixture<KineticsThermoCurveFit> & _chem;

    //! Primary thermo object.
    std::unique_ptr<Thermo> _thermo;

    //! For some Thermo types, we also need to cache a NASAEvaluator.
    std::unique_ptr<Antioch::NASAEvaluator<libMesh::Real,KineticsThermoCurveFit> > _nasa_evaluator;

    std::unique_ptr<AntiochKinetics<KineticsThermoCurveFit> > _kinetics;

    // Temperature should be clipped to positive values for
    // stability's sake
    libMesh::Real _clipped_T;

    // We'll clip temperature at a user-specified minimum if
    // requested.  Some reaction equations give us NaNs at 0K too.
    const libMesh::Real _minimum_T;

    std::unique_ptr<Antioch::TempCache<libMesh::Real> > _temp_cache;

    //! Helper method for managing _temp_cache
    /*! T *MUST* be pass-by-reference because of the structure
      of Antioch::TempCache! */
    void check_and_reset_temp_cache( const libMesh::Real& T );

  private:

    // StatMechThermodynamics methods
    libMesh::Real specialized_cp( const libMesh::Real & T,
                                  const libMesh::Real /*P*/,
                                  const std::vector<libMesh::Real> & Y,
                                  const Antioch::StatMechThermodynamics<libMesh::Real> & thermo )
    {
      return thermo.cp( T, T, Y );
    }

    void specialized_cp_s( const libMesh::Real & T,
                           const libMesh::Real /*P*/,
                           const std::vector<libMesh::Real> & Y,
                           const Antioch::StatMechThermodynamics<libMesh::Real> & thermo,
                           std::vector<libMesh::Real> & cp_s )
    {
      libmesh_assert_equal_to(_chem.n_species(), Y.size());
      libmesh_assert_equal_to(Y.size(), cp_s.size());

      for(unsigned int s = 0; s < Y.size(); s++)
        cp_s[s] = thermo.cv(s, T, T) + this->R(s);
    }

    libMesh::Real specialized_cv( const libMesh::Real & T,
                                  const libMesh::Real /*P*/,
                                  const std::vector<libMesh::Real> & Y,
                                  const Antioch::StatMechThermodynamics<libMesh::Real> & thermo )
    { return thermo.cv( T, T, Y ); }

    libMesh::Real specialized_h_s( const libMesh::Real & T, unsigned int species,
                                   const Antioch::StatMechThermodynamics<libMesh::Real> & thermo )
    { return thermo.h_tot( species, T ) + this->_chem.h_stat_mech_ref_correction(species); }



    // IdealGasThermo methods
    libMesh::Real specialized_cp( const libMesh::Real & T,
                                  const libMesh::Real /*P*/,
                                  const std::vector<libMesh::Real> & Y,
                                  const Antioch::IdealGasThermo<KineticsThermoCurveFit,libMesh::Real> & /*thermo*/ )
    {
      this->check_and_reset_temp_cache(T);
      return this->_nasa_evaluator->cp( *_temp_cache, Y );
    }

    void specialized_cp_s( const libMesh::Real & T,
                           const libMesh::Real /*P*/,
                           const std::vector<libMesh::Real> & Y,
                           const Antioch::IdealGasThermo<KineticsThermoCurveFit,libMesh::Real> & /*thermo*/,
                           std::vector<libMesh::Real> & cp_s )
    {
      libmesh_assert_equal_to(_chem.n_species(), Y.size());
      libmesh_assert_equal_to(Y.size(), cp_s.size());

      this->check_and_reset_temp_cache(T);
      for(unsigned int s = 0; s < Y.size(); s++)
        cp_s[s] = _nasa_evaluator->cp( *_temp_cache, s );
    }

    libMesh::Real specialized_cv( const libMesh::Real & T,
                                  const libMesh::Real /*P*/,
                                  const std::vector<libMesh::Real> & Y,
                                  const Antioch::IdealGasThermo<KineticsThermoCurveFit,libMesh::Real> & /*thermo*/ )
    {
      this->check_and_reset_temp_cache(T);
      return this->_nasa_evaluator->cv( *_temp_cache, Y );
    }

    libMesh::Real specialized_h_s( const libMesh::Real & T, unsigned int species,
                                   const Antioch::IdealGasThermo<KineticsThermoCurveFit,libMesh::Real> & /*thermo*/ )
    {
      this->check_and_reset_temp_cache(T);
      return this->_nasa_evaluator->h( *_temp_cache, species );
    }

  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename KineticsThermoCurveFit, typename Thermo>
  inline
  libMesh::Real
  AntiochEvaluator<KineticsThermoCurveFit,Thermo>::M( unsigned int species ) const
  {
    return _chem.M(species);
  }

  template<typename KineticsThermoCurveFit, typename Thermo>
  inline
  libMesh::Real
  AntiochEvaluator<KineticsThermoCurveFit,Thermo>::
  M_mix( const std::vector<libMesh::Real>& mass_fractions ) const
  {
    return _chem.M_mix(mass_fractions);
  }

  template<typename KineticsThermoCurveFit, typename Thermo>
  inline
  libMesh::Real
  AntiochEvaluator<KineticsThermoCurveFit,Thermo>::R( unsigned int species ) const
  {
    return _chem.R(species);
  }

  template<typename KineticsThermoCurveFit, typename Thermo>
  inline
  libMesh::Real
  AntiochEvaluator<KineticsThermoCurveFit,Thermo>::
  R_mix( const std::vector<libMesh::Real>& mass_fractions ) const
  {
    return _chem.R_mix(mass_fractions);
  }

  template<typename KineticsThermoCurveFit, typename Thermo>
  inline
  libMesh::Real AntiochEvaluator<KineticsThermoCurveFit,Thermo>::X( unsigned int species,
                                                                    libMesh::Real M,
                                                                    libMesh::Real mass_fraction ) const
  {
    return _chem.X(species,M,mass_fraction);
  }

  template<typename KineticsThermoCurveFit, typename Thermo>
  inline
  void AntiochEvaluator<KineticsThermoCurveFit,Thermo>::X( libMesh::Real M,
                                                           const std::vector<libMesh::Real>& mass_fractions,
                                                           std::vector<libMesh::Real>& mole_fractions ) const
  {
    _chem.X(M,mass_fractions,mole_fractions);
  }

  template<typename KineticsThermoCurveFit, typename Thermo>
  inline
  unsigned int
  AntiochEvaluator<KineticsThermoCurveFit,Thermo>::species_index( const std::string& species_name ) const
  {
    return _chem.species_index(species_name);
  }

  template<typename KineticsThermoCurveFit, typename Thermo>
  inline
  std::string
  AntiochEvaluator<KineticsThermoCurveFit,Thermo>::species_name( unsigned int species_index ) const
  {
    return _chem.species_name(species_index);
  }

  template<typename KineticsThermoCurveFit, typename Thermo>
  inline
  void AntiochEvaluator<KineticsThermoCurveFit,Thermo>::check_and_reset_temp_cache( const libMesh::Real& T )
  {
    _clipped_T = std::max(_minimum_T, T);

    // We can't compare T because it's a reference, so we may have already
    // changed it upstream. So, we compare the next cheapest thing.
    if( _temp_cache->T2 != _clipped_T*_clipped_T )
      _temp_cache.reset( new Antioch::TempCache<libMesh::Real>(_clipped_T) );
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_EVALUATOR_H
