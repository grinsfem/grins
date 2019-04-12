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


#ifndef GRINS_PARSED_PRESSURE_H
#define GRINS_PARSED_PRESSURE_H

//GRINS
#include "grins/common.h"
#include "grins/parsed_property_base.h"
#include "grins/parameter_user.h"

// libMesh
#include "libmesh/getpot.h"
class GetPot;

namespace GRINS
{
  //! Class to manage spatially varying pressure parameter from parsed function input
  class ParsedPressure : public ParsedPropertyBase<ParsedPressure>,
                         public ParameterUser
  {
  public:

    //! This will parse the input for <section>/pressure
    ParsedPressure( const GetPot & input, const std::string & section );

    ParsedPressure() = delete;

    virtual ~ParsedPressure() = default;
  };

  inline
  ParsedPressure::ParsedPressure( const GetPot & input, const std::string & section )
    : ParsedPropertyBase(),
      ParameterUser("ParsedPressure")
  {
    std::string var = section+"/pressure";

    if( !input.have_variable(var) )
      libmesh_error_msg("Error: Must supply input pressure using key "+var+"!\n");

    this->set_parameter(this->_func, input, var, "DIE!" );
  }

} // end namespace GRINS


#endif // GRINS_PARSED_PRESSURE_H
