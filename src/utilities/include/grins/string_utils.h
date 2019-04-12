//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_STRING_UTILS_H
#define GRINS_STRING_UTILS_H

// libMesh
#include "libmesh/libmesh_common.h"

// C++
#include <sstream>
#include <string>
#include <vector>

namespace GRINS
{
  namespace StringUtilities
  {
    template <typename T>
    inline
    T string_to_T(const std::string& input)
    {
      std::istringstream converter(input);
      T returnval;
      converter >> returnval;
      if (converter.fail())
        libmesh_error();
      return returnval;
    }

    template <typename T>
    inline
    std::string T_to_string(const T input)
    {
      std::stringstream converter;
      converter << input;
      if (converter.fail())
        libmesh_error();
      return converter.str();
    }

    /*!
      Split on colon, and return name, int value pair.
      Taken from FIN-S for XML parsing.
    */
    inline
    std::pair<std::string, int> split_string_int_on_colon(const std::string &token)
    {
      std::pair<std::string, int> ret = std::make_pair(std::string(), 0);
      std::string::size_type colon_position = token.find(":");
      libmesh_assert (colon_position != std::string::npos);
      ret.first  = token.substr(0, colon_position);
      ret.second = std::atoi(token.substr(colon_position + 1).c_str());
      return ret;
    }


    /*!
      Split on colon, and return name, double value pair.
      Taken from FIN-S for XML parsing.
    */
    inline
    std::pair<std::string, double> split_string_double_on_colon(const std::string &token)
    {
      std::pair<std::string, double> ret = std::make_pair(std::string(), 0.0);
      std::string::size_type colon_position = token.find(":");
      libmesh_assert (colon_position != std::string::npos);
      ret.first  = token.substr(0, colon_position);
      ret.second = std::atof(token.substr(colon_position + 1).c_str());
      return ret;
    }

    void split_string( const std::string& input,
                       const std::string& delimiter,
                       std::vector<std::string>& results );

    void split_string_real( const std::string& input,
                            const std::string& delimiter,
                            std::vector<libMesh::Real>& results );
  } // end namespace StringUtilities
} // end namespace GRINS

#endif // GRINS_STRING_UTILS_H
