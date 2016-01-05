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


#ifndef GRINS_THERMO_PRESSURE_VARIABLE_H
#define GRINS_THERMO_PRESSURE_VARIABLE_H

// GRINS
#include "grins/variables_base.h"

// libMesh forward declarations
class GetPot;
namespace libMesh
{
  class FEMSystem;
}

namespace GRINS
{
  class ThermoPressureVariable : public VariablesBase
  {
  public:

    ThermoPressureVariable( const GetPot& input );
    ~ThermoPressureVariable(){};

    virtual void init( libMesh::FEMSystem* system );

    VariableIndex p0_var() const;

  private:

    ThermoPressureVariable();

  };

  inline
  VariableIndex ThermoPressureVariable::p0_var() const
  {
    return _vars[0];
  }

} // end namespace GRINS

#endif // GRINS_THERMO_PRESSURE_VARIABLE_H
