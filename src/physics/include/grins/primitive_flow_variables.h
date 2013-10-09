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


#ifndef GRINS_PRIMITIVE_FLOW_VARIABLES_H
#define GRINS_PRIMITIVE_FLOW_VARIABLES_H

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

  class PrimitiveFlowVariables
  {
  public:

    PrimitiveFlowVariables( const GetPot& input );
    ~PrimitiveFlowVariables();

    virtual void init( libMesh::FEMSystem* system );

    VariableIndex u_var() const;
    VariableIndex v_var() const;
    VariableIndex w_var() const;
    VariableIndex p_var() const;

  protected:

    //! Indices for each (owned) variable;
    VariableIndex _u_var; /* Index for x-velocity field */
    VariableIndex _v_var; /* Index for y-velocity field */
    VariableIndex _w_var; /* Index for z-velocity field */
    VariableIndex _p_var; /* Index for pressure field */

    //! Names of each (owned) variable in the system
    std::string _u_var_name, _v_var_name, _w_var_name, _p_var_name;

  private:

    PrimitiveFlowVariables();

  };

  inline
  VariableIndex PrimitiveFlowVariables::u_var() const
  {
    return _u_var;
  }

  inline
  VariableIndex PrimitiveFlowVariables::v_var() const
  {
    return _v_var;
  }

  inline
  VariableIndex PrimitiveFlowVariables::w_var() const
  {
    return _w_var;
  }

  inline
  VariableIndex PrimitiveFlowVariables::p_var() const
  {
    return _p_var;
  }

} // end namespace GRINS

#endif //GRINS_PRIMITIVE_FLOW_VARIABLES_H
