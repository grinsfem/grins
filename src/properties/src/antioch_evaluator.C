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


#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// C++
#include <limits>

// This class
#include "grins/antioch_evaluator.h"

// GRINS
#include "grins/antioch_mixture.h"
#include "grins/cached_values.h"
#include "grins/antioch_mixture_builder_base.h"

namespace GRINS
{
  template<typename KineticsThermoCurveFit, typename Thermo>
  AntiochEvaluator<KineticsThermoCurveFit,Thermo>::
  AntiochEvaluator( const AntiochMixture<KineticsThermoCurveFit>& mixture )
    : _chem( mixture ),
      _nasa_evaluator( new Antioch::NASAEvaluator<libMesh::Real,KineticsThermoCurveFit>(mixture.nasa_mixture()) ),
      _kinetics( new AntiochKinetics<KineticsThermoCurveFit>(mixture) ),
      _minimum_T( mixture.minimum_T() ),
      _temp_cache( new Antioch::TempCache<libMesh::Real>(1.0) )
  {
    AntiochMixtureBuilderBase builder;
    _thermo = builder.build_gas_thermo<KineticsThermoCurveFit,Thermo>
      ( mixture.chemical_mixture(), mixture.nasa_mixture() );
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

} // end namespace GRINS

#endif //GRINS_HAVE_ANTIOCH
