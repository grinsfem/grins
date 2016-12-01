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

// This class
#include "grins/builder_helper.h"

// GRINS
#include "grins/string_utils.h"
#include "grins/variables_parsing.h"

// libMesh
#include "libmesh/getpot.h"

// C++
#include <algorithm>

namespace GRINS
{
  void BuilderHelper::parse_var_sections( const GetPot& input,
                                          std::set<std::string>& sections )
  {
    std::vector<std::string> sec_vec;
    parse_var_sections_vector(input,sec_vec);

    // Now convert populated vector to a set
    for(std::vector<std::string>::const_iterator it = sec_vec.begin();
        it != sec_vec.end(); ++it )
      sections.insert(*it);
  }

  void BuilderHelper::parse_var_sections_vector( const GetPot& input,
                                                 std::vector<std::string>& sections )
  {
    if( !input.have_section(VariablesParsing::variables_section()) )
       libmesh_error_msg("ERROR: Could not find "+VariablesParsing::variables_section()+" section!");

    // We need to extract all the Variable subsections from the input file
    sections = input.get_subsection_names(VariablesParsing::variables_section());

    // Make sure we found some variable subsections
    if( sections.empty() )
       libmesh_error_msg("ERROR: Did not find any Variable subsections!");
  }
} // end namespace GRINS
