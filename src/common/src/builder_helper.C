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
#include "grins/builder_helper.h"

// GRINS
#include "grins/string_utils.h"
#include "grins/variables_parsing.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  void BuilderHelper::parse_var_sections( const GetPot& input,
                                          std::set<std::string>& sections ) const
  {
    if( !input.have_section(VariablesParsing::variables_section()) )
       libmesh_error_msg("ERROR: Could not find "+VariablesParsing::variables_section()+" section!");

    // We need to extract all the Variable sections from the input file
    // We'll populate the relevant sections in var_sections
    /*! \todo This would probably be a good function to add to libMesh::GetPot */
    std::vector<std::string> all_sections = input.get_section_names();
    for( std::vector<std::string>::const_iterator s = all_sections.begin();
         s < all_sections.end(); ++s )
      {
        // First check that it contains "Variable" as the first slot
        if( s->find(VariablesParsing::variables_section()) == 0 )
          {
            // Now check it only has 2 elements when we split on "/"
            std::vector<std::string> split_str;
            StringUtilities::split_string(*s, "/", split_str );

            // Our Variable should be the second part of the split
            if( split_str.size() == 2 )
              {
                // Make sure we don't already have that section
                if( sections.find(split_str[1]) != sections.end() )
                  libmesh_error_msg("ERROR: Found duplicate Variable section "+split_str[1]+"!");

                sections.insert( split_str[1] );
              }
          }
      }

    // Make sure we found some variable subsections
    if( sections.empty() )
       libmesh_error_msg("ERROR: Did not find any Variable subsections!");
  }
} // end namespace GRINS
