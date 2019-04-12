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

#ifndef GRINS_SIMULATION_PARSING_H
#define GRINS_SIMULATION_PARSING_H

// C++
#include <string>

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  class SimulationParsing
  {
  public:

    SimulationParsing(){};

    ~SimulationParsing(){};

    static bool have_restart( const GetPot& input )
    { return input.have_variable( SimulationParsing::restart_input_option() ); }

    static std::string restart_file( const GetPot& input )
    { return input( SimulationParsing::restart_input_option(), "none" ); }

  private:

    static std::string restart_input_option()
    { return "restart-options/restart_file"; }
  };

} // end namespace GRINS

#endif // GRINS_SIMULATION_PARSING_H
