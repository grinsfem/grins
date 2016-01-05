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

// This class
#include "grins/species_mass_fracs_variables.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

// GRINS
#include "grins/variable_name_defaults.h"
#include "grins/materials_parsing.h"

namespace GRINS
{
  SpeciesMassFractionsVariables::SpeciesMassFractionsVariables( const GetPot& input,
                                                                const std::string& material_name )
  {
    MaterialsParsing::parse_species_varnames(input, material_name, _species_var_names);
  }

  void SpeciesMassFractionsVariables::init( libMesh::FEMSystem* system )
  {
    unsigned int n_species = this->n_species();

#ifndef NDEBUG
    for( unsigned int s = 0; s < n_species; s++)
      libmesh_assert( system->has_variable( _species_var_names[s] ) );
#endif

    _species_vars.resize(n_species);
    for( unsigned int s = 0; s < n_species; s++)
      _species_vars[s] = system->variable_number( _species_var_names[s] );
  }

} // end namespace GRINS
