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

// This class
#include "grins/antioch_mixture_builder_base.h"

// GRINS
#include "grins/common.h"
#include "grins/materials_parsing.h"

// Antioch
#include "antioch/default_filename.h"
#include "antioch/read_reaction_set_data.h"
#include "antioch/xml_parser.h"

namespace GRINS
{
  std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> >
  AntiochMixtureBuilderBase::build_chem_mix( const GetPot & input, const std::string & material )
  {
    std::vector<std::string> species_list;
    this->build_species_names(input,material,species_list);

    std::string prefix(this->antioch_prefix(material));

    bool verbose_antioch_read = input(prefix+"/verbose_read",false);

    std::string species_data_filename = input(prefix+"/species_data", "default" );
    if( species_data_filename == std::string("default") )
      species_data_filename = Antioch::DefaultInstallFilename::chemical_mixture();

    std::string vibration_data_filename = input(prefix+"/vibration_data", "default" );
    if( vibration_data_filename == std::string("default") )
      vibration_data_filename = Antioch::DefaultInstallFilename::vibrational_data();

    std::string electronic_data_filename = input(prefix+"/electronic_data", "default" );
    if( electronic_data_filename == std::string("default") )
      electronic_data_filename = Antioch::DefaultInstallFilename::electronic_data();

    // By default, Antioch is using its ASCII parser. We haven't added more options yet.
    return std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> >
      ( new Antioch::ChemicalMixture<libMesh::Real>( species_list,
                                                     verbose_antioch_read,
                                                     species_data_filename,
                                                     vibration_data_filename,
                                                     electronic_data_filename ) );
  }

  std::unique_ptr<Antioch::ReactionSet<libMesh::Real> >
  AntiochMixtureBuilderBase::build_reaction_set( const GetPot & input, const std::string & material,
                                                 const Antioch::ChemicalMixture<libMesh::Real> & chem_mix )
  {
    std::unique_ptr<Antioch::ReactionSet<libMesh::Real> >
      reaction_set( new Antioch::ReactionSet<libMesh::Real>(chem_mix) );

    bool verbose_read = input("screen-options/verbose_kinetics_read", false );

    std::string prefix(this->antioch_prefix(material));
    std::string gas_mixture_option(prefix+"/gas_mixture");

    Antioch::ParsingType parsing_type = this->get_antioch_parsing_type(input,material);

    switch(parsing_type)
      {
        // We're actually going to use Antioch's XML parser to parse the kinetics,
        // but we distinguish ASCII because we need a separate GRINS input option
        // to get the kinetics file. In the XML case, everything is in one XML file
        // whereas with the ASCII, we have \infty files...
      case(Antioch::ASCII):
        {
          std::string old_option("Materials/"+material+"/GasMixture/kinetics_data");
          std::string new_option(prefix+"/kinetics_data");

          std::string filename;
          if( input.have_variable(old_option) )
            {
              // Deprecated message is in this function call
              filename = MaterialsParsing::parse_chemical_kinetics_datafile_name( input, material );
            }
          else if( input.have_variable(new_option) )
            filename = input(new_option,std::string("DIE!"));

          else
            libmesh_error_msg("ERROR: Could not valid input for "+new_option+"!\n");

          if( input.have_variable(gas_mixture_option) )
            {
              std::string gas_mixture = input(gas_mixture_option,"DIE!");
              Antioch::XMLParser<libMesh::Real> parser(filename, gas_mixture, verbose_read);
              Antioch::read_reaction_set_data<libMesh::Real>( verbose_read, *reaction_set, &parser );
            }
          else
            {
              std::string msg = "WARNING: Option "+gas_mixture_option+" not found!\n";
              msg += "         The first kinetics data set in the file "+filename+"\n";
              msg += "         will be used. This is DEPRECATED behavior. In the future,\n";
              msg += "         the option "+gas_mixture_option+" will be required.\n";
              grins_warning(msg);

              Antioch::read_reaction_set_data<libMesh::Real>( filename, verbose_read, *reaction_set, Antioch::XML );
            }
          break;
        }
      case(Antioch::XML):
        {
          std::unique_ptr<Antioch::XMLParser<libMesh::Real> > parser = this->build_xml_parser(input,material);
          Antioch::read_reaction_set_data<libMesh::Real>( verbose_read, *reaction_set, parser.get() );
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



    return reaction_set;
  }

  void AntiochMixtureBuilderBase::build_species_names( const GetPot & input, const std::string & material,
                                                       std::vector<std::string> & species_names)
  {
    // Clear out anything the user might've put in there.
    species_names.clear();

    Antioch::ParsingType parsing_type = this->get_antioch_parsing_type(input,material);

    std::string prefix(this->antioch_prefix(material));

    switch(parsing_type)
      {
      case(Antioch::ASCII):
        {
          std::string old_option("Materials/"+material+"/GasMixture/species");
          std::string new_option(prefix+"/species");

          if( input.have_variable(old_option) )
            this->parse_chemical_species(input,old_option,species_names);

          else if( input.have_variable(new_option) )
            this->parse_chemical_species(input,new_option,species_names);

          else
            {
              std::string msg = "ERROR: Could not find valid entry for "+new_option+" !\n";
              msg += "       Using the ASCII parser with Antioch, you must explicitly\n";
              msg += "       specify the list of species using the input option\n";
              msg += "       "+new_option+" .\n";
              msg += "       We strongly recommend switching to XML format.\n";

              libmesh_error_msg(msg);
            }
          break;
        }
      case(Antioch::XML):
        {
          std::unique_ptr<Antioch::XMLParser<libMesh::Real> > parser = this->build_xml_parser(input,material);
          species_names = parser->species_list();
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

  void AntiochMixtureBuilderBase::parse_chemical_species( const GetPot & input,
                                                          const std::string & option,
                                                          std::vector<std::string>& species_names )
  {

    MaterialsParsing::check_for_input_option(input,option);

    // Read variable naming info
    unsigned int n_species = input.vector_variable_size(option);

    species_names.reserve(n_species);
    for( unsigned int i = 0; i < n_species; i++ )
      species_names.push_back( input( option, "DIE!", i ) );

  }

  Antioch::ParsingType
  AntiochMixtureBuilderBase::get_antioch_parsing_type( const GetPot & input, const std::string & material ) const
  {
    Antioch::ParsingType parsing_type = Antioch::ASCII;

    std::string prefix(this->antioch_prefix(material));
    std::string chem_data_option(prefix+"/"+MaterialsParsing::chemical_data_option());

    // Check for old_option. If this is specified, then we are doing ASCII.
    // We only really want to support XML and ChemKin parsing interfaces.
    std::string old_option("Materials/"+material+"/GasMixture/species");
    if( input.have_variable(old_option) )
      {
        parsing_type = Antioch::ASCII;
        std::string warning = "WARNING: Specifying "+old_option+" is DEPRECATED!\n";
        warning += "         Instead use the option "+prefix+"/species\n";
        warning += "         or, preferred, use XML or ChemKin formats so that species\n";
        warning += "         are parsed automatically.\n";
        grins_warning(warning);
      }

    // If not the above, then we're new style. Check the suffix of the chemical data file.
    // If .dat --> ASCII
    // If .xml --> XML
    // If .chemkin --> CHEMKIN
    else if( input.have_variable(prefix+"/"+MaterialsParsing::chemical_data_option()) )
      {
        std::string chem_data_filename = input(chem_data_option, std::string("DIE!") );
        if( chem_data_filename.find(".dat") != chem_data_filename.npos )
          parsing_type = Antioch::ASCII;
        else if( chem_data_filename.find(".xml") != chem_data_filename.npos )
          parsing_type = Antioch::XML;
        else if( chem_data_filename.find(".chemkin") != chem_data_filename.npos )
          parsing_type = Antioch::CHEMKIN;
        else
          {
            std::string msg = "ERROR: Could not find valid input specification for ";
            msg += prefix+"/"+MaterialsParsing::chemical_data_option()+"\n";
            msg += "       Expected a filename with .dat, .xml, or .chemkin suffix.\n";
            libmesh_error_msg(msg);
          }
      }
    else
      {
        std::string msg = "ERROR: Could not find valid input specification for ";
        msg += prefix+"/"+MaterialsParsing::chemical_data_option()+"\n";
        libmesh_error_msg(msg);
      }

    return parsing_type;
  }

  std::unique_ptr<Antioch::XMLParser<libMesh::Real> >
  AntiochMixtureBuilderBase::build_xml_parser( const GetPot & input, const std::string & material ) const
  {
    std::string prefix(this->antioch_prefix(material));

    std::string chem_data_option(prefix+"/"+MaterialsParsing::chemical_data_option());
    std::string chem_data_filename = input(chem_data_option,std::string("DIE!"));

    std::string gas_mixture_option(prefix+"/gas_mixture");
    MaterialsParsing::check_for_input_option(input,gas_mixture_option);

    std::string gas_mixture = input(gas_mixture_option,"DIE!");

    return std::unique_ptr<Antioch::XMLParser<libMesh::Real> >
      ( new Antioch::XMLParser<libMesh::Real>(chem_data_filename, gas_mixture, false) );
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
