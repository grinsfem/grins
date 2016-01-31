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

#ifndef GRINS_TURBULENCE_VARIABLES_H
#define GRINS_TURBULENCE_VARIABLES_H

// libMesh forward declarations
class GetPot;

// GRINS
#include "grins/single_variable.h"

namespace GRINS
{
  class TurbulenceVariables : public SingleVariable
  {
  public:

    TurbulenceVariables( const GetPot& input )
      : SingleVariable(input,
                       this->old_var_name(),
                       this->subsection(),
                       this->default_name())
    {}

    ~TurbulenceVariables(){};

    VariableIndex nu() const;

  private:

    std::string old_var_name() const
    { return "turbulent_viscosity"; }

    std::string subsection() const
    { return "TurbulentViscosity"; }

    std::string default_name() const
    { return "nu"; }

    TurbulenceVariables();

  };

  inline
  VariableIndex TurbulenceVariables::nu() const
  {
    return _vars[0];
  }

} // end namespace GRINS

#endif // GRINS_TURBULENCE_VARIABLES_H
