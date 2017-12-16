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

// This class
#include "grins/variables_parsing.h"

namespace GRINS
{
  std::string VariablesParsing::section_parse_var_name( const GetPot & input,
                                                        const std::string & input_subsection_name,
                                                        const std::string & input_var_name,
                                                        const std::string & default_var_name,
                                                        const SECTION_TYPE section_type )
  {
    std::string header;
    switch(section_type)
      {
      case(PHYSICS):
        {
          header = "Physics";
          break;
        }
      case(QOI):
        {
          header = "QoI";
          break;
        }
      default:
        libmesh_error_msg("ERROR: Invalid section value!");
      }

    return input(header+"/"+input_subsection_name+"/"+input_var_name, default_var_name );
  }
} // end namespace GRINS
