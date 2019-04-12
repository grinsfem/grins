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


#ifndef GRINS_PARSED_CONDUCTIVITY_H
#define GRINS_PARSED_CONDUCTIVITY_H

//GRINS
#include "grins/parameter_user.h"
#include "grins/parsed_property_base.h"

class GetPot;

namespace GRINS
{
  class ParsedConductivity : public ParsedPropertyBase<ParsedConductivity>,
                             public ParameterUser
  {
  public:

    //! Constructor with specified material
    /*! Will look in the input file for [Materials/material/ThermalConductivity/value]
      for the value of viscosity. */
    ParsedConductivity( const GetPot& input, const std::string& material );

    //! Deprecated constructor
    ParsedConductivity( const GetPot& input );
    ~ParsedConductivity();

  private:

    ParsedConductivity();

  };
} // end namespace GRINS

#endif // GRINS_PARSED_CONDUCTIVITY_H
