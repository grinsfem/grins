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


#ifndef GRINS_WILKE_ANTIOCH_TRANSPORT_MIXTURE_H
#define GRINS_WILKE_ANTIOCH_TRANSPORT_MIXTURE_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// GRINS
#include "grins/antioch_mixture.h"
#include "grins/property_types.h"
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/getpot.h"

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/vector_utils.h"
#include "antioch/cea_curve_fit.h"
#include "antioch/nasa_evaluator.h"
#include "antioch/ideal_gas_thermo.h"
#include "antioch/stat_mech_thermo.h"
#include "antioch/mixture_averaged_transport_mixture.h"
#include "antioch/mixture_viscosity.h"
#include "antioch/mixture_conductivity.h"
#include "antioch/mixture_diffusion.h"
#include "antioch/sutherland_viscosity.h"
#include "antioch/blottner_viscosity.h"
#include "antioch/eucken_thermal_conductivity.h"
#include "antioch/constant_lewis_diffusivity.h"

#include "antioch/sutherland_parsing.h"
#include "antioch/blottner_parsing.h"
#include "antioch/eucken_thermal_conductivity_building.h"
#include "antioch/constant_lewis_diffusivity_building.h"

#ifdef ANTIOCH_HAVE_GSL

#include "antioch/kinetics_theory_viscosity.h"
#include "antioch/kinetics_theory_thermal_conductivity.h"
#include "antioch/molecular_binary_diffusion.h"
#include "antioch/gsl_spliner.h"
#include "antioch/kinetics_theory_viscosity_building.h"
#include "antioch/kinetics_theory_thermal_conductivity_building.h"

#endif // ANTIOCH_HAVE_GSL

namespace GRINS
{
  //! Wrapper class for storing state for computing Wilke transport properties using Antioch
  /*!
    This class is expected to be constructed *before* threads have been forked and will
    live during the whole program.
    By default, Antioch is working in SI units. Note that this documentation will always
    be built regardless if Antioch is included in the GRINS build or not. Check configure
    output to confirm that Antioch was included in the build.
  */
  template<typename KineticsThermoCurveFit, typename Thermo, typename Viscosity, typename Conductivity, typename Diffusivity>
  class AntiochMixtureAveragedTransportMixture : public AntiochMixture<KineticsThermoCurveFit>
  {
  public:

    //! Deprecated Constructor
    AntiochMixtureAveragedTransportMixture( const GetPot& input, const std::string& material );

    //! Constructor with user-built objects
    AntiochMixtureAveragedTransportMixture( std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> > & chem_mixture,
                                            std::unique_ptr<Antioch::ReactionSet<libMesh::Real> > & reaction_set,
                                            std::unique_ptr<Antioch::NASAThermoMixture<libMesh::Real,KineticsThermoCurveFit> > & kinetics_thermo_mix,
                                            std::unique_ptr<Thermo> & gas_thermo,
                                            std::unique_ptr<Antioch::TransportMixture<libMesh::Real> > & trans_mix,
                                            std::unique_ptr<Antioch::MixtureAveragedTransportMixture<libMesh::Real> > & wilke_mix,
                                            std::unique_ptr<Antioch::MixtureViscosity<Viscosity,libMesh::Real> > & visc,
                                            std::unique_ptr<Antioch::MixtureConductivity<Conductivity,libMesh::Real> > & conductivity,
                                            std::unique_ptr<Antioch::MixtureDiffusion<Diffusivity,libMesh::Real> > & diff,
                                            libMesh::Real min_T = -std::numeric_limits<libMesh::Real>::max(),
                                            bool clip_negative_rho = false );

    virtual ~AntiochMixtureAveragedTransportMixture(){}

    const Antioch::MixtureAveragedTransportMixture<libMesh::Real>& wilke_mixture() const;

    const Antioch::MixtureViscosity<Viscosity,libMesh::Real>& viscosity() const;

    const Antioch::MixtureConductivity<Conductivity,libMesh::Real>& conductivity() const;

    const Antioch::MixtureDiffusion<Diffusivity,libMesh::Real>& diffusivity() const;

    typedef AntiochChemistry ChemistryParent;

  protected:

    std::unique_ptr<Antioch::TransportMixture<libMesh::Real> > _trans_mixture;

    std::unique_ptr<Antioch::MixtureAveragedTransportMixture<libMesh::Real> > _wilke_mixture;

    std::unique_ptr<Thermo> _thermo;

    std::unique_ptr<Antioch::MixtureViscosity<Viscosity,libMesh::Real> > _viscosity;

    std::unique_ptr<Antioch::MixtureConductivity<Conductivity,libMesh::Real> > _conductivity;

    std::unique_ptr<Antioch::MixtureDiffusion<Diffusivity,libMesh::Real> > _diffusivity;

  private:

    AntiochMixtureAveragedTransportMixture();

  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename KT, typename T, typename V, typename C, typename D>
  inline
  const Antioch::MixtureAveragedTransportMixture<libMesh::Real>&
  AntiochMixtureAveragedTransportMixture<KT,T,V,C,D>::wilke_mixture() const
  {
    return *(_wilke_mixture.get());
  }

  template<typename KT, typename T, typename V, typename C, typename D>
  inline
  const Antioch::MixtureViscosity<V,libMesh::Real>&
  AntiochMixtureAveragedTransportMixture<KT,T,V,C,D>::viscosity() const
  {
    return *_viscosity.get();
  }

  template<typename KT, typename T, typename V, typename C, typename D>
  inline
  const Antioch::MixtureConductivity<C,libMesh::Real>&
  AntiochMixtureAveragedTransportMixture<KT,T,V,C,D>::conductivity() const
  {
    return *_conductivity.get();
  }

  template<typename KT, typename T, typename V, typename C, typename D>
  inline
  const Antioch::MixtureDiffusion<D,libMesh::Real>&
  AntiochMixtureAveragedTransportMixture<KT,T,V,C,D>::diffusivity() const
  {
    return *_diffusivity.get();
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_WILKE_ANTIOCH_TRANSPORT_MIXTURE_H
