//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_SOLID_MECHANICS_VARIABLES_H
#define GRINS_SOLID_MECHANICS_VARIABLES_H

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
  class SolidMechanicsVariables
  {
  public:

    SolidMechanicsVariables( const GetPot& input );
    virtual ~SolidMechanicsVariables();

    virtual void init( libMesh::FEMSystem* system );

    VariableIndex u_var() const;
    VariableIndex v_var() const;
    VariableIndex w_var() const;

  protected:

    VariableIndex _u_var;
    VariableIndex _v_var;
    VariableIndex _w_var;

    std::string _u_var_name, _v_var_name, _w_var_name;

  };

} // end namespace GRINS

#endif // GRINS_SOLID_MECHANICS_VARIABLES_H
