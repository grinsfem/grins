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

#ifndef GRINS_SOLVER_PARSING_H
#define GRINS_SOLVER_PARSING_H

// C++
#include <string>

// Forward declarations
class GetPot;

namespace GRINS
{
  class SolverParsing
  {
  public:

    SolverParsing(){};

    ~SolverParsing(){};

    static std::string solver_type(const GetPot& input);

    static bool is_transient( const GetPot& input );

    static void dup_solver_option_check( const GetPot& input,
                                         const std::string& option1,
                                         const std::string& option2 );
  };

} // end namespace GRINS

#endif // GRINS_SOLVER_PARSING_H
