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


#ifndef GRINS_ANTIOCH_MIXTURE_BUILDER_BASE_H
#define GRINS_ANTIOCH_MIXTURE_BUILDER_BASE_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// GRINS
#include "grins/materials_parsing.h"
#include "grins/antioch_mixture.h"

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/vector_utils.h"
#include "antioch/chemical_mixture.h"
#include "antioch/reaction_set.h"
#include "antioch/nasa_mixture.h"
#include "antioch/cea_curve_fit.h"
#include "antioch/nasa_mixture_parsing.h"

// libMesh
#include "libmesh/auto_ptr.h" // libMesh::UniquePtr
#include "libmesh/getpot.h"

// C++
#include <string>

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  //! Base class building Antioch mixture wrappers
  /*! This class only worries about building the kinetics
      and the thermo associated with kinetics. Subclasses
      will handle thermo and transport. */
  class AntiochMixtureBuilderBase
  {
  public:
    AntiochMixtureBuilderBase(){}
    ~AntiochMixtureBuilderBase(){}

    libMesh::UniquePtr<Antioch::ChemicalMixture<libMesh::Real> >
    build_chem_mix( const GetPot & input, const std::string & material );

    libMesh::UniquePtr<Antioch::ReactionSet<libMesh::Real> >
    build_reaction_set( const GetPot & input, const std::string & material,
                        const Antioch::ChemicalMixture<libMesh::Real> & chem_mix );

    template<typename KineticsThermoCurveFit>
    libMesh::UniquePtr<Antioch::NASAThermoMixture<libMesh::Real,KineticsThermoCurveFit> >
    build_nasa_thermo_mix( const GetPot & input, const std::string & material,
                           const Antioch::ChemicalMixture<libMesh::Real> & chem_mix );


    template<typename KineticsThermoCurveFit>
    libMesh::UniquePtr<AntiochMixture<KineticsThermoCurveFit> >
    build_antioch_mixture(const GetPot & input, const std::string & material );

    libMesh::Real parse_min_T( const GetPot & input, const std::string & material )
    {
      return input( "Materials/"+material+"/GasMixture/Antioch/minimum_T",
                    -std::numeric_limits<libMesh::Real>::max() );
    }

    bool parse_clip_negative_rho( const GetPot & input, const std::string & material )
    {
      return input( "Materials/"+material+"/GasMixture/Antioch/clip_negative_rho", false);
    }

  protected:

    void parse_nasa_data
    ( Antioch::NASAThermoMixture<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> > & nasa_mixture,
      const GetPot & input, const std::string & material)
    {
      std::string cea_data_filename = input( "Materials/"+material+"/GasMixture/Antioch/cea_data", "default" );

      if( cea_data_filename == std::string("default") )
        cea_data_filename = Antioch::DefaultInstallFilename::thermo_data();

      Antioch::read_nasa_mixture_data( nasa_mixture, cea_data_filename, Antioch::ASCII, true );
    }

  };

  template<typename KineticsThermoCurveFit>
  inline
  libMesh::UniquePtr<Antioch::NASAThermoMixture<libMesh::Real,KineticsThermoCurveFit> >
  AntiochMixtureBuilderBase::build_nasa_thermo_mix( const GetPot & input, const std::string & material,
                                                    const Antioch::ChemicalMixture<libMesh::Real> & chem_mix )
  {
    libMesh::UniquePtr<Antioch::NASAThermoMixture<libMesh::Real,KineticsThermoCurveFit> >
      nasa_mixture( new Antioch::NASAThermoMixture<libMesh::Real,KineticsThermoCurveFit>(chem_mix) );

    this->parse_nasa_data( *nasa_mixture, input, material);

    return nasa_mixture;
  }

  template<typename KineticsThermoCurveFit>
  inline
  libMesh::UniquePtr<AntiochMixture<KineticsThermoCurveFit> >
  AntiochMixtureBuilderBase::build_antioch_mixture(const GetPot & input, const std::string & material )
  {
    libMesh::UniquePtr<Antioch::ChemicalMixture<libMesh::Real> > chem_mixture =
      this->build_chem_mix(input,material);

    libMesh::UniquePtr<Antioch::ReactionSet<libMesh::Real> > reaction_set =
      this->build_reaction_set(input,material,*chem_mixture);

    libMesh::UniquePtr<Antioch::NASAThermoMixture<libMesh::Real,KineticsThermoCurveFit> > nasa_mixture =
      this->build_nasa_thermo_mix<KineticsThermoCurveFit>(input,material,*chem_mixture);

    libMesh::Real min_T = this->parse_min_T(input,material);
    bool clip_negative_rho = this->parse_clip_negative_rho(input,material);

    return libMesh::UniquePtr<AntiochMixture<KineticsThermoCurveFit> >
      ( new AntiochMixture<KineticsThermoCurveFit>
        (chem_mixture,reaction_set,nasa_mixture,min_T,clip_negative_rho) );
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_MIXTURE_BUILDER_BASE_H
