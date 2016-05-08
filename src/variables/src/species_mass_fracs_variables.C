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
#include "libmesh/fem_system.h"

namespace GRINS
{
  void SpeciesMassFractionsVariables::vars_from_species_prefix( libMesh::FEMSystem* system )
  {
    std::vector<unsigned int> all_var_nums;
    system->get_all_variable_numbers(all_var_nums);

    for( std::vector<unsigned int>::const_iterator it = all_var_nums.begin();
         it < all_var_nums.end(); ++it )
      {
        const std::string& var_name = system->variable_name(*it);

        if( var_name.find(_prefix) != std::string::npos )
          {
            _var_names.push_back(var_name);
            _vars.push_back(*it);
          }
      }
  }

} // end namespace GRINS
