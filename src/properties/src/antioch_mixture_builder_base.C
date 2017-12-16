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

// This class
#include "grins/antioch_mixture_builder_base.h"

// GRINS
#include "grins/materials_parsing.h"

// Antioch
#include "antioch/default_filename.h"
#include "antioch/read_reaction_set_data.h"

namespace GRINS
{
  std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> >
  AntiochMixtureBuilderBase::build_chem_mix( const GetPot & input, const std::string & material )
  {
    std::vector<std::string> species_list;
    MaterialsParsing::parse_chemical_species(input,material,species_list);

    bool verbose_antioch_read = input("Materials/"+material+"/GasMixture/Antioch/verbose_read",false);

    std::string species_data_filename = input("Materials/"+material+"/GasMixture/Antioch/species_data", "default" );
    if( species_data_filename == std::string("default") )
      species_data_filename = Antioch::DefaultInstallFilename::chemical_mixture();

    std::string vibration_data_filename = input("Materials/"+material+"/GasMixture/Antioch/vibration_data", "default" );
    if( vibration_data_filename == std::string("default") )
      vibration_data_filename = Antioch::DefaultInstallFilename::vibrational_data();

    std::string electronic_data_filename = input("Materials/"+material+"/GasMixture/Antioch/electronic_data", "default" );
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

    std::string kinetics_data_filename = MaterialsParsing::parse_chemical_kinetics_datafile_name( input, material );

    bool verbose_read = input("screen-options/verbose_kinetics_read", false );

    Antioch::read_reaction_set_data_xml<libMesh::Real>( kinetics_data_filename, verbose_read, *reaction_set );

    return reaction_set;
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
