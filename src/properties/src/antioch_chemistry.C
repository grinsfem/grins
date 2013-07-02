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
#include "grins/antioch_chemistry.h"

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

    _antioch_gas.reset( new Antioch::ChemicalMixture<libMesh::Real>( species_list ) );

    return;
  }

  AntiochChemistry::~AntiochChemistry()
  {
    return;
  }
  
  std::string AntiochChemistry::species_name( unsigned int species_index ) const
  {
    libmesh_assert_less(species_index, _antioch_gas->n_species());

    std::string name = "dummy";

    for( std::map<std::string,unsigned int>::const_iterator it = _antioch_gas->active_species_name_map().begin();
         it != _antioch_gas->active_species_name_map().end(); it++ )
      {
        if( it->second == species_index )
          {
            name = it->first;
          }
      }

    if( name == std::string("dummy") )
      {
        std::cerr << "Error: Could not find a species name for the given index!"
                  << std::endl;
        libmesh_error();
      }

    return name;
  }

}// end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
