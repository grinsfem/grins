//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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
#include "grins/antioch_chemistry.h"

// Antioch
#include "antioch/default_filename.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  AntiochChemistry::AntiochChemistry( const GetPot& input )
    : _antioch_gas(NULL)
  {
    if( !input.have_variable("Physics/Chemistry/species") )
      {
        std::cerr << "Error: Must specify species list to use Antioch." << std::endl;
        libmesh_error();
      }

    unsigned int n_species = input.vector_variable_size("Physics/Chemistry/species");
    std::vector<std::string> species_list(n_species);

    for( unsigned int s = 0; s < n_species; s++ )
      {
        species_list[s] = input( "Physics/Chemistry/species", "DIE!", s );
      }

    bool verbose_antioch_read = input("Physics/Antioch/verbose_read",false);

    std::string species_data_filename = input("Physics/Antioch/species_data", "default" );
    if( species_data_filename == std::string("default") )
      species_data_filename = Antioch::DefaultInstallFilename::chemical_mixture();

    std::string vibration_data_filename = input("Physics/Antioch/vibration_data", "default" );
    if( vibration_data_filename == std::string("default") )
      vibration_data_filename = Antioch::DefaultInstallFilename::vibrational_data();

    std::string electronic_data_filename = input("Physics/Antioch/electronic_data", "default" );
    if( electronic_data_filename == std::string("default") )
      electronic_data_filename = Antioch::DefaultInstallFilename::electronic_data();

    // By default, Antioch is using its ASCII parser. We haven't added more options yet.
    _antioch_gas.reset( new Antioch::ChemicalMixture<libMesh::Real>( species_list,
                                                                     verbose_antioch_read,
                                                                     species_data_filename,
                                                                     vibration_data_filename,
                                                                     electronic_data_filename ) );

    return;
  }

  AntiochChemistry::~AntiochChemistry()
  {
    return;
  }
  
  std::string AntiochChemistry::species_name( unsigned int species_index ) const
  {
    libmesh_assert_less(species_index, _antioch_gas->n_species());

    return _antioch_gas->species_inverse_name_map().find(species_index)->second;
  }

}// end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
