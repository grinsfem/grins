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
#include "antioch/xml_parser.h"

// libMesh
#include "libmesh/auto_ptr.h" // std::unique_ptr
#include "libmesh/getpot.h"

// C++
#include <string>

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  enum ThermoEnum {NASA7 = 0, NASA9, CEA, INVALID};

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

    //! Returns a unique_ptr to a fully constructed Antioch::ChemicalMixture
    std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> >
    build_chem_mix( const GetPot & input, const std::string & material );

    //! Returns a unique_ptr to a fully constructed Antioch::ReactionSet
    std::unique_ptr<Antioch::ReactionSet<libMesh::Real> >
    build_reaction_set( const GetPot & input, const std::string & material,
                        const Antioch::ChemicalMixture<libMesh::Real> & chem_mix );

    ThermoEnum get_thermo_type( const GetPot & input, const std::string & material ) const;

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
    ( const GetPot & input, const std::string & material, Antioch::ParsingType parsing_type,
      Antioch::NASAThermoMixture<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> > & nasa_mixture )
    {
      switch(parsing_type)
      {
      case(Antioch::ASCII):
        {
          std::string prefix(this->antioch_prefix(material));
          std::string cea_data_filename = input( prefix+"/cea_data", "default" );

          if( cea_data_filename == std::string("default") )
            cea_data_filename = Antioch::DefaultInstallFilename::thermo_data();

          Antioch::read_nasa_mixture_data( nasa_mixture, cea_data_filename, Antioch::ASCII, false );

          break;
        }
      case(Antioch::XML):
        {
          libmesh_error_msg("ERROR: XML Parsing of CEA not implemented. Use NASA7 or NASA9 instead!");
          break;
        }
      case(Antioch::CHEMKIN):
        {
          libmesh_error_msg("ERROR: ChemKin Parsing of CEA not implemented. Use NASA7!");
          break;
        }
      default:
        libmesh_error_msg("ERROR: Invalid Antioch parsing type!");
      }
    }

    void parse_nasa_data
    ( const GetPot & input, const std::string & material, Antioch::ParsingType parsing_type,
      Antioch::NASAThermoMixture<libMesh::Real,Antioch::NASA7CurveFit<libMesh::Real> > & nasa_mixture )
    {
      switch(parsing_type)
      {
      case(Antioch::ASCII):
        {
          libmesh_error_msg("ERROR: ASCII Parsing of NASA7 not implemented. Use NASA9 instead for ASCII parser!");
          break;
        }
      case(Antioch::XML):
        {
          std::unique_ptr<Antioch::XMLParser<libMesh::Real> > parser = this->build_xml_parser(input,material);
          libmesh_assert(parser->is_nasa7_curve_fit_type());
          parser->read_thermodynamic_data(nasa_mixture);
          break;
        }
      case(Antioch::CHEMKIN):
        {
          libmesh_not_implemented();
          break;
        }
      default:
        libmesh_error_msg("ERROR: Invalid Antioch parsing type!");
      }
    }

    void parse_nasa_data
    ( const GetPot & input, const std::string & material, Antioch::ParsingType parsing_type,
      Antioch::NASAThermoMixture<libMesh::Real,Antioch::NASA9CurveFit<libMesh::Real> > & nasa_mixture )
    {
      switch(parsing_type)
      {
      case(Antioch::ASCII):
        {
          libmesh_error_msg("ERROR: ASCII Parsing of NASA7 not implemented. Use NASA9 instead for ASCII parser!");
          break;
        }
      case(Antioch::XML):
        {
          std::unique_ptr<Antioch::XMLParser<libMesh::Real> > parser = this->build_xml_parser(input,material);
          libmesh_assert(!parser->is_nasa7_curve_fit_type());
          parser->read_thermodynamic_data(nasa_mixture);
          break;
        }
      case(Antioch::CHEMKIN):
        {
          libmesh_not_implemented();
          break;
        }
      default:
        libmesh_error_msg("ERROR: Invalid Antioch parsing type!");
      }
    }

    //! Helper function for parsing the chemical species
    /*! The user-provided vector will populated with the chemical
        species names in the input file. "Material/"+material+"/GasMixture/species"
        is the option name. */
    void parse_chemical_species( const GetPot & input,
                                 const std::string & option,
                                 std::vector<std::string>& species_names );

    //! Determine the input file type the user is using for Antioch
    /*! There are three valid options: ASCII, XML, and CHEMKIN. We will use a convention on the
        suffix of the chemical_data file name to ascertain the parsing type. Once we know the
        parsing type from this function, other functions can act accordingly. */
    Antioch::ParsingType get_antioch_parsing_type( const GetPot & input, const std::string & material ) const;

    //! Helper function to encapsulate parsing prefix
    std::string antioch_prefix( const std::string & material ) const
    { return std::string("Materials/"+material+"/GasMixture/Antioch"); }

    //! Helper function to build an XML parser when we need it.
    /*! Builds an XMLParser using the chemical_data input option and the gas_mixture option.
      Currently, not every *single* XMLParser is initialized this way, but most are.
      Hopefully in the future we can stream line this a bit, but for now with
      the current interface, we're building an XML parser several times, so let's
      put that in one place. */
    std::unique_ptr<Antioch::XMLParser<libMesh::Real> >
    build_xml_parser( const GetPot & input, const std::string & material ) const;

  };

  template<typename NASACurveFit>
  inline
  std::unique_ptr<Antioch::NASAThermoMixture<libMesh::Real,NASACurveFit> >
  AntiochMixtureBuilderBase::build_nasa_thermo_mix( const GetPot & input, const std::string & material,
                                                    const Antioch::ChemicalMixture<libMesh::Real> & chem_mix )
  {
    std::unique_ptr<Antioch::NASAThermoMixture<libMesh::Real,NASACurveFit> >
      nasa_mixture( new Antioch::NASAThermoMixture<libMesh::Real,NASACurveFit>(chem_mix) );

    Antioch::ParsingType parsing_type = this->get_antioch_parsing_type(input,material);

    this->parse_nasa_data( input, material, parsing_type, *nasa_mixture);

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

  inline
  ThermoEnum AntiochMixtureBuilderBase::get_thermo_type( const GetPot & input,
                                                         const std::string & material ) const
  {
    ThermoEnum thermo_type = ThermoEnum::INVALID;

    Antioch::ParsingType parsing_type = this->get_antioch_parsing_type(input,material);
    switch(parsing_type)
      {
      case(Antioch::ASCII):
        {
          // ASCII Parsing is only CEA, so it'd better be CEA
          thermo_type = ThermoEnum::CEA;
          break;
        }
        case(Antioch::XML):
        {
          // XML can be NASA7 or NASA9, so we query the XML file to see which is in there.
          std::unique_ptr<Antioch::XMLParser<libMesh::Real> > parser = this->build_xml_parser(input,material);
          if( parser->is_nasa7_curve_fit_type() )
            thermo_type = ThermoEnum::NASA7;
          else
            thermo_type = ThermoEnum::NASA9;
          break;
        }
      case(Antioch::CHEMKIN):
        {
          // ChemKin is NASA7 by definition of ChemKin format.
          thermo_type = ThermoEnum::NASA7;
          break;
        }
      default:
        libmesh_error_msg("ERROR: Invalid Antioch parsing type!");
      }

    return thermo_type;
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_MIXTURE_BUILDER_BASE_H
