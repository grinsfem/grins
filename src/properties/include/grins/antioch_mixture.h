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
#include "antioch/nasa_mixture.h"
#include "antioch/cea_curve_fit.h"
#include "antioch/reaction_set.h"

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  //! Wrapper class for storing state for Antioch thermo and kinetics
  /*!
    This class handles caching needed state for Antioch kinetics, and the thermodynamics
    required for the kinetics evaluation. Currently, we only support NASA type curve fits
    for required thermodynamic evaluations, but we template on the curve fit type.

    This class is expected to be constructed *before* threads have been forked and will
    live during the whole program.
    By default, Antioch is working in SI units. Note that this documentation will always
    be built regardless if Antioch is included in the GRINS build or not. Check configure
    output to confirm that Antioch was included in the build.
  */
  template <typename KineticsThermoCurveFit>
  class AntiochMixture : public AntiochChemistry
  {
  public:

    //! Deprecated Constructor
    AntiochMixture( const GetPot& input, const std::string& material );


    //! Constructor with user-built objects
    /*! This constructor expects the user to pass in a ChemicalMixture, ReactionSet, and NASAThermoMixture
      object already built; this class will take ownership of the pointer. */
    AntiochMixture( std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> > & chem_mixture,
                    std::unique_ptr<Antioch::ReactionSet<libMesh::Real> > & reaction_set,
                    std::unique_ptr<Antioch::NASAThermoMixture<libMesh::Real,KineticsThermoCurveFit> > & nasa_mixture,
                    libMesh::Real min_T = -std::numeric_limits<libMesh::Real>::max(),
                    bool clip_negative_rho = false );

    virtual ~AntiochMixture(){};

    // Registers all parameters in this physics and in its property
    // classes
    virtual void register_parameter
    ( const std::string & param_name,
      libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
      const;

    const Antioch::ReactionSet<libMesh::Real>& reaction_set() const;

    const Antioch::NASAThermoMixture<libMesh::Real,KineticsThermoCurveFit> & nasa_mixture() const;

    libMesh::Real h_stat_mech_ref_correction( unsigned int species ) const;

    // Returns the minimum temperature at which reactions will be
    // evaluated, in Kelvin
    libMesh::Real minimum_T() const;

    // Returns true iff negative species densities are clipped to zero
    // when calculating reaction rates.
    bool clip_negative_rho() const;

  protected:

    std::unique_ptr<Antioch::ReactionSet<libMesh::Real> > _reaction_set;

    std::unique_ptr<Antioch::NASAThermoMixture<libMesh::Real,KineticsThermoCurveFit> > _nasa_mixture;

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
  template <typename KineticsThermoCurveFit>
  inline
  const Antioch::ReactionSet<libMesh::Real>&
  AntiochMixture<KineticsThermoCurveFit>::reaction_set() const
  {
    return *_reaction_set.get();
  }

  template <typename KineticsThermoCurveFit>
  inline
  const Antioch::NASAThermoMixture<libMesh::Real,KineticsThermoCurveFit> &
  AntiochMixture<KineticsThermoCurveFit>::nasa_mixture() const
  {
    return *_nasa_mixture.get();
  }

  template <typename KineticsThermoCurveFit>
  inline
  libMesh::Real
  AntiochMixture<KineticsThermoCurveFit>::h_stat_mech_ref_correction( unsigned int species ) const
  {
    return _h_stat_mech_ref_correction[species];
  }

  template <typename KineticsThermoCurveFit>
  inline
  libMesh::Real AntiochMixture<KineticsThermoCurveFit>::minimum_T() const
  {
    return _minimum_T;
  }

  template <typename KineticsThermoCurveFit>
  inline
  bool AntiochMixture<KineticsThermoCurveFit>::clip_negative_rho() const
  {
    return _clip_negative_rho;
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_MIXTURE_H
