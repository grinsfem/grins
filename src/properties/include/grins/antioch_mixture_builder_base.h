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


#ifndef GRINS_ANTIOCH_MIXTURE_BUILDER_BASE_H
#define GRINS_ANTIOCH_MIXTURE_BUILDER_BASE_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// GRINS
#include "grins/materials_parsing.h"
#include "grins/property_types.h"
#include "grins/antioch_mixture.h"

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/vector_utils.h"
#include "antioch/chemical_mixture.h"
#include "antioch/reaction_set.h"
#include "antioch/nasa_mixture.h"
#include "antioch/cea_curve_fit.h"
#include "antioch/nasa_mixture_parsing.h"
#include "antioch/ideal_gas_thermo.h"
#include "antioch/stat_mech_thermo.h"

// libMesh
#include "libmesh/auto_ptr.h" // std::unique_ptr
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

    //! Populates the passed in vector of strings with all the species names
    /*! This method handles the details of how we parse the species names using Antioch's
        parsing methods. The vector of strings will be sized accordingly, the user
        does not need to worry about correct sizing. */
    void
    build_species_names( const GetPot & input, const std::string & material,
                         std::vector<std::string> & species_names);

    std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> >
    build_chem_mix( const GetPot & input, const std::string & material );

    std::unique_ptr<Antioch::ReactionSet<libMesh::Real> >
    build_reaction_set( const GetPot & input, const std::string & material,
                        const Antioch::ChemicalMixture<libMesh::Real> & chem_mix );

    template<typename NASACurveFit>
    std::unique_ptr<Antioch::NASAThermoMixture<libMesh::Real,NASACurveFit> >
    build_nasa_thermo_mix( const GetPot & input, const std::string & material,
                           const Antioch::ChemicalMixture<libMesh::Real> & chem_mix );

    template<typename NASACurveFit, typename Thermo>
    std::unique_ptr<Thermo>
    build_gas_thermo( const Antioch::ChemicalMixture<libMesh::Real> & chem_mix,
                      const Antioch::NASAThermoMixture<libMesh::Real,NASACurveFit> & nasa_mix )
    { return specialized_build_gas_thermo<NASACurveFit>( chem_mix, nasa_mix, thermo_type<Thermo>() ); }

    template<typename NASACurveFit>
    std::unique_ptr<AntiochMixture<NASACurveFit> >
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

    template<typename NASACurveFit>
    std::unique_ptr<Antioch::StatMechThermodynamics<libMesh::Real> >
    specialized_build_gas_thermo( const Antioch::ChemicalMixture<libMesh::Real> & chem_mix,
                                  const Antioch::NASAThermoMixture<libMesh::Real,NASACurveFit> & /*nasa_mix*/,
                                  thermo_type<Antioch::StatMechThermodynamics<libMesh::Real> > )
    {
      return std::unique_ptr<Antioch::StatMechThermodynamics<libMesh::Real> >
        (new Antioch::StatMechThermodynamics<libMesh::Real>(chem_mix));
    }

    template<typename NASACurveFit>
    std::unique_ptr<Antioch::IdealGasThermo<NASACurveFit,libMesh::Real> >
    specialized_build_gas_thermo( const Antioch::ChemicalMixture<libMesh::Real> & chem_mix,
                                  const Antioch::NASAThermoMixture<libMesh::Real,NASACurveFit> & nasa_mix,
                                  thermo_type<Antioch::IdealGasThermo<NASACurveFit,libMesh::Real> > )
    {
      return std::unique_ptr<Antioch::IdealGasThermo<NASACurveFit,libMesh::Real> >
        (new Antioch::IdealGasThermo<NASACurveFit,libMesh::Real>(nasa_mix,chem_mix));
    }

    void parse_nasa_data
    ( Antioch::NASAThermoMixture<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> > & nasa_mixture,
      const GetPot & input, const std::string & material)
    {
      std::string cea_data_filename = input( "Materials/"+material+"/GasMixture/Antioch/cea_data", "default" );

      if( cea_data_filename == std::string("default") )
        cea_data_filename = Antioch::DefaultInstallFilename::thermo_data();

      Antioch::read_nasa_mixture_data( nasa_mixture, cea_data_filename, Antioch::ASCII, true );
    }

    void parse_nasa_data
    ( Antioch::NASAThermoMixture<libMesh::Real,Antioch::NASA7CurveFit<libMesh::Real> > & /*nasa_mixture*/,
      const GetPot & /*input*/, const std::string & /*material*/)
    {
      libmesh_not_implemented();
    }

    void parse_nasa_data
    ( Antioch::NASAThermoMixture<libMesh::Real,Antioch::NASA9CurveFit<libMesh::Real> > & /*nasa_mixture*/,
      const GetPot & /*input*/, const std::string & /*material*/)
    {
      libmesh_not_implemented();
    }

    //! Helper function for parsing the chemical species
    /*! The user-provided vector will populated with the chemical
        species names in the input file. "Material/"+material+"/GasMixture/species"
        is the option name. */
    void parse_chemical_species( const GetPot & input,
                                 const std::string & material,
                                 std::vector<std::string>& species_names );

  };

  template<typename NASACurveFit>
  inline
  std::unique_ptr<Antioch::NASAThermoMixture<libMesh::Real,NASACurveFit> >
  AntiochMixtureBuilderBase::build_nasa_thermo_mix( const GetPot & input, const std::string & material,
                                                    const Antioch::ChemicalMixture<libMesh::Real> & chem_mix )
  {
    std::unique_ptr<Antioch::NASAThermoMixture<libMesh::Real,NASACurveFit> >
      nasa_mixture( new Antioch::NASAThermoMixture<libMesh::Real,NASACurveFit>(chem_mix) );

    this->parse_nasa_data( *nasa_mixture, input, material);

    return nasa_mixture;
  }

  template<typename NASACurveFit>
  inline
  std::unique_ptr<AntiochMixture<NASACurveFit> >
  AntiochMixtureBuilderBase::build_antioch_mixture(const GetPot & input, const std::string & material )
  {
    std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> > chem_mixture =
      this->build_chem_mix(input,material);

    std::unique_ptr<Antioch::ReactionSet<libMesh::Real> > reaction_set =
      this->build_reaction_set(input,material,*chem_mixture);

    std::unique_ptr<Antioch::NASAThermoMixture<libMesh::Real,NASACurveFit> > nasa_mixture =
      this->build_nasa_thermo_mix<NASACurveFit>(input,material,*chem_mixture);

    libMesh::Real min_T = this->parse_min_T(input,material);
    bool clip_negative_rho = this->parse_clip_negative_rho(input,material);

    return std::unique_ptr<AntiochMixture<NASACurveFit> >
      ( new AntiochMixture<NASACurveFit>
        (chem_mixture,reaction_set,nasa_mixture,min_T,clip_negative_rho) );
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_MIXTURE_BUILDER_BASE_H
