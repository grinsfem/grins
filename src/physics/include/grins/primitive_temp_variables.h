//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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

#ifndef GRINS_PRIMITIVE_TEMP_VARIABLES_H
#define GRINS_PRIMITIVE_TEMP_VARIABLES_H

// libMesh forward declarations
class GetPot;
namespace libMesh
{
  class FEMSystem;
}

// GRINS
#include "grins/var_typedefs.h"

namespace GRINS
{
  class PrimitiveTempVariables
  {
  public:

    PrimitiveTempVariables( const GetPot& input );
    ~PrimitiveTempVariables();

    virtual void init( libMesh::FEMSystem* system );

    VariableIndex T_var() const;

  protected:

    //! Indices for each variable;
    VariableIndex _T_var; /* Index for temperature field */

    //! Names of each variable in the system
    std::string _T_var_name;

  private:

    PrimitiveTempVariables();

  };

  inline
  VariableIndex PrimitiveTempVariables::T_var() const
  {
    return _T_var;
  }

} // end namespace GRINS

#endif // GRINS_PRIMITIVE_TEMP_VARIABLES_H
