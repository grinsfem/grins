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

#ifndef GRINS_BUILDER_HELPER_H
#define GRINS_BUILDER_HELPER_H

// C++
#include <set>
#include <string>
#include <vector>

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  //! This class contains common functions used by various Builders throughout GRINS
  class BuilderHelper
  {
  public:
    BuilderHelper(){}
    ~BuilderHelper(){}

    //! Parses the input file for [Variables] first-level subsections
    /*! The expected format for the Variables is
        \code{.unparsed}
        [Variables]
           [./Displacement]
              fe_family = '...'
              order = '...'
              names = '...'
           [../]
           [./Velocity]
              ...
           [../]
        [../]
        \endcode
        For the example above, this function will fill 'sections' with "Displacement"
        and "Velocity".

        The result is return in an std::set so the actual ordering from the input file
        is lost.
    */
    static
    void parse_var_sections( const GetPot& input,
                             std::set<std::string>& sections );

    //! The same as parse_var_sections, except the result is returned in an std::vector.
    /*! This allows the maintaining of the order of the Variable subsection in the input
        file. This is important for things like setting up the Variables in the System
        so that the user can control the order the variables are added to the system.*/
    static
    void parse_var_sections_vector( const GetPot& input,
                                    std::vector<std::string>& sections );
  };
} // end namespace GRINS

#endif // GRINS_BUILDER_HELPER_H
