//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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

// This class
#include "grins/antioch_kinetics.h"

// GRINS
#include "grins/antioch_mixture.h"

// Antioch
#include "antioch/temp_cache.h"
#include "antioch/vector_utils.h"

namespace GRINS
{
  AntiochKinetics::AntiochKinetics( const AntiochMixture& mixture )
    : _antioch_mixture( mixture ),
      _antioch_kinetics( mixture.reaction_set(), 0 ),
      _antioch_cea_thermo( mixture.cea_mixture() )
  {}

  void AntiochKinetics::omega_dot( const libMesh::Real& T,
                                   const libMesh::Real rho,
                                   const std::vector<libMesh::Real>& mass_fractions,
                                   std::vector<libMesh::Real>& omega_dot )
  {
    Antioch::TempCache<libMesh::Real> temp_cache(T);
    this->omega_dot(temp_cache,rho,mass_fractions,omega_dot);
  }

  void AntiochKinetics::omega_dot( const Antioch::TempCache<libMesh::Real>& temp_cache,
                                   const libMesh::Real rho,
                                   const std::vector<libMesh::Real>& mass_fractions,
                                   std::vector<libMesh::Real>& omega_dot )
  {
    const unsigned int n_species = _antioch_mixture.n_species();

    libmesh_assert_equal_to( mass_fractions.size(), n_species );
    libmesh_assert_equal_to( omega_dot.size(), n_species );

    std::vector<libMesh::Real> h_RT_minus_s_R(n_species, 0.0);
    std::vector<libMesh::Real> molar_densities(n_species, 0.0);

    _antioch_cea_thermo.h_RT_minus_s_R( temp_cache, h_RT_minus_s_R );

    _antioch_mixture.molar_densities( rho, mass_fractions, molar_densities );

    bool have_density = false;
    for (unsigned int i=0; i != n_species; ++i)
      if (molar_densities[i] <= 0)
        molar_densities[i] = 0;
      else
        have_density = true;

    if (have_density)
      _antioch_kinetics.compute_mass_sources( temp_cache.T,
                                              molar_densities,
                                              h_RT_minus_s_R,
                                              omega_dot );
    else
      std::fill(omega_dot.begin(), omega_dot.end(), 0.0);
  }

}// end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
