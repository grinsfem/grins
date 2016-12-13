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


#ifndef GRINS_ANTIOCH_MIXTURE_H
#define GRINS_ANTIOCH_MIXTURE_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// GRINS
#include "grins/antioch_chemistry.h"
#include "grins/property_types.h"

// libMesh
#include "libmesh/libmesh_common.h"

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/vector_utils.h"
#include "antioch/chemical_mixture.h"
#include "antioch/cea_mixture.h"
#include "antioch/reaction_set.h"

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  //! Wrapper class for storing state for Antioch thermo and kinetics
  /*!
    This class is expected to be constructed *before* threads have been forked and will
    live during the whole program.
    By default, Antioch is working in SI units. Note that this documentation will always
    be built regardless if Antioch is included in the GRINS build or not. Check configure
    output to confirm that Antioch was included in the build.
   */
  class AntiochMixture : public AntiochChemistry
  {
  public:

    AntiochMixture( const GetPot& input, const std::string& material );

    virtual ~AntiochMixture(){};

    // Registers all parameters in this physics and in its property
    // classes
    virtual void register_parameter
      ( const std::string & param_name,
        libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
    const;

    const Antioch::ReactionSet<libMesh::Real>& reaction_set() const;

    const Antioch::CEAThermoMixture<libMesh::Real>& cea_mixture() const;

    libMesh::Real h_stat_mech_ref_correction( unsigned int species ) const;

    // Returns the minimum temperature at which reactions will be
    // evaluated, in Kelvin
    libMesh::Real minimum_T() const;

    // Returns true iff negative species densities are clipped to zero
    // when calculating reaction rates.
    bool clip_negative_rho() const;

  protected:

    libMesh::UniquePtr<Antioch::ReactionSet<libMesh::Real> > _reaction_set;

    libMesh::UniquePtr<Antioch::CEAThermoMixture<libMesh::Real> > _cea_mixture;

    std::vector<libMesh::Real> _h_stat_mech_ref_correction;

    void build_stat_mech_ref_correction();

    // Users can specify a minimum temperature at which to evaluate
    // reaction rates.  Some reaction equations give us NaNs at 0
    // Kelvin or less, and solution "ringing" can result in those
    // unphysical temperatures for some formulations near strong
    // fronts.
    //
    // By default, temperatures are not clipped to any minimum.
    libMesh::Real _minimum_T;

    // Users can request that negative densities be clipped to zero
    // when evaluating reaction rates.  Some reaction equations fail
    // badly with negative densities input, and solution "ringing" can
    // result in those unphysical densities for some formulations near
    // strong fronts.
    //
    // By default, negative densities are not clipped to zero.
    bool _clip_negative_rho;

  private:

    AntiochMixture();

  };

  /* ------------------------- Inline Functions -------------------------*/
  inline
  const Antioch::ReactionSet<libMesh::Real>& AntiochMixture::reaction_set() const
  {
    return *_reaction_set.get();
  }

  inline
  const Antioch::CEAThermoMixture<libMesh::Real>& AntiochMixture::cea_mixture() const
  {
    return *_cea_mixture.get();
  }

  inline
  libMesh::Real AntiochMixture::h_stat_mech_ref_correction( unsigned int species ) const
  {
    return _h_stat_mech_ref_correction[species];
  }

  inline
  libMesh::Real AntiochMixture::minimum_T() const
  {
    return _minimum_T;
  }

  inline
  bool AntiochMixture::clip_negative_rho() const
  {
    return _clip_negative_rho;
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_MIXTURE_H
