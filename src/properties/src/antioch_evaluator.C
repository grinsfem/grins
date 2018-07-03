//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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


#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// C++
#include <limits>

// This class
#include "grins/antioch_evaluator.h"

// GRINS
#include "grins/antioch_mixture.h"
#include "grins/cached_values.h"

namespace GRINS
{
  template<typename KineticsThermoCurveFit, typename Thermo>
  AntiochEvaluator<KineticsThermoCurveFit,Thermo>::
  AntiochEvaluator( const AntiochMixture<KineticsThermoCurveFit>& mixture )
    : _chem( mixture ),
      _kinetics( new AntiochKinetics<KineticsThermoCurveFit>(mixture) ),
      _minimum_T( mixture.minimum_T() ),
      _temp_cache( new Antioch::TempCache<libMesh::Real>(1.0) )
  {
    this->build_thermo( mixture );
  }

  template<typename KineticsThermoCurveFit, typename Thermo>
  void AntiochEvaluator<KineticsThermoCurveFit,Thermo>::
  omega_dot( const libMesh::Real& T, libMesh::Real rho,
             const std::vector<libMesh::Real> mass_fractions,
             std::vector<libMesh::Real>& omega_dot )
  {
    this->check_and_reset_temp_cache(T);

    _kinetics->omega_dot( *(_temp_cache.get()), rho, mass_fractions, omega_dot );
  }

  template<>
  libMesh::Real
  AntiochEvaluator<Antioch::CEACurveFit<libMesh::Real>,Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real> >::
  cp( const libMesh::Real& T,
      const libMesh::Real /*P*/,
      const std::vector<libMesh::Real>& Y )
  {
    this->check_and_reset_temp_cache(T);
    return this->_nasa_evaluator->cp( *(_temp_cache.get()), Y );
  }

  template<>
  void AntiochEvaluator<Antioch::CEACurveFit<libMesh::Real>,Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real> >::
  cp_s( const libMesh::Real& T,
	const libMesh::Real /*P*/,
	const std::vector<libMesh::Real>& Y,
        std::vector<libMesh::Real>& Cp_s )
  {
    libmesh_assert_equal_to(Y.size(), Cp_s.size());
    this->check_and_reset_temp_cache(T);
    for(unsigned int species = 0;species < Y.size();species++)
      Cp_s[species] = _nasa_evaluator->cp( *(_temp_cache.get()), species );
    return;
  }

  template<>
  libMesh::Real AntiochEvaluator<Antioch::CEACurveFit<libMesh::Real>,Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real> >::
  cv( const libMesh::Real& T,
      const libMesh::Real /*P*/,
      const std::vector<libMesh::Real>& Y )
  {
    this->check_and_reset_temp_cache(T);
    return this->_nasa_evaluator->cv( *(_temp_cache.get()), Y );
  }

  template<>
  libMesh::Real AntiochEvaluator<Antioch::CEACurveFit<libMesh::Real>,Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real> >::
  h_s( const libMesh::Real& T, unsigned int species )
  {
    this->check_and_reset_temp_cache(T);

    return this->_nasa_evaluator->h( *(_temp_cache.get()), species );;
  }

  template<>
  libMesh::Real AntiochEvaluator<Antioch::CEACurveFit<libMesh::Real>,Antioch::StatMechThermodynamics<libMesh::Real> >::cp( const libMesh::Real& T,
                                                                                                                           const libMesh::Real /*P*/,
                                                                                                                           const std::vector<libMesh::Real>& Y )
  {
    return this->_thermo->cp( T, T, Y );
  }

  template<>
  void AntiochEvaluator<Antioch::CEACurveFit<libMesh::Real>,Antioch::StatMechThermodynamics<libMesh::Real> >::cp_s( const libMesh::Real& T,
														    const libMesh::Real /*P*/,
														    const std::vector<libMesh::Real>& Y,
														    std::vector<libMesh::Real>& Cp_s)
  {
    //will need to figure out before pushing
    for(unsigned int species = 0;species< Y.size();species++)
      Cp_s[species] = _thermo->cv(species, T, T) + this->R(species);
    return;
  }


  template<>
  libMesh::Real AntiochEvaluator<Antioch::CEACurveFit<libMesh::Real>,Antioch::StatMechThermodynamics<libMesh::Real> >::cv( const libMesh::Real& T,
                                                                                                                           const libMesh::Real /*P*/,
                                                                                                                           const std::vector<libMesh::Real>& Y )
  {
    return this->_thermo->cv( T, T, Y );
  }

  template<>
  libMesh::Real AntiochEvaluator<Antioch::CEACurveFit<libMesh::Real>,Antioch::StatMechThermodynamics<libMesh::Real> >::h_s( const libMesh::Real& T, unsigned int species )
  {
    return this->_thermo->h_tot( species, T ) + _chem.h_stat_mech_ref_correction(species);
  }

} // end namespace GRINS

#endif //GRINS_HAVE_ANTIOCH
