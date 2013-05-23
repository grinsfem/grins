//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// This class
#include "grins/antioch_mixture.h"

// libMesh
#include "libmesh/getpot.h"

// Antioch
#include "antioch/read_reaction_set_data_xml.h"
#include "antioch/cea_mixture_ascii_parsing.h"

namespace GRINS
{
  AntiochMixture::AntiochMixture( const GetPot& input )
    : _antioch_gas(NULL),
      _reaction_set(NULL),
      _cea_mixture(NULL)
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

    _antioch_gas.reset( new Antioch::ChemicalMixture<libMesh::Real>( species_list ) );

    _reaction_set.reset( new Antioch::ReactionSet<libMesh::Real>( (*_antioch_gas.get()) ) );

    if( !input.have_variable("Physics/Chemistry/chem_file") )
      {
        std::cerr << "Error: Must specify XML chemistry file to use Antioch." << std::endl;
        libmesh_error();
      }

    std::string xml_filename = input( "Physics/Chemistry/chem_file", "DIE!");
    bool verbose_read = input("screen-options/verbose_kinetics_read", false );
      
    Antioch::read_reaction_set_data_xml<libMesh::Real>( xml_filename, verbose_read, *_reaction_set.get() );

    _cea_mixture.reset( new Antioch::CEAThermoMixture<libMesh::Real>( (*_antioch_gas.get()) ) );

    Antioch::read_cea_mixture_data_ascii_default( *_cea_mixture.get() );

    return;
  }

  AntiochMixture::~AntiochMixture()
  {
    return;
  }

}// end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
