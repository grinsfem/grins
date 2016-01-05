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


#ifndef GRINS_PRIMITIVE_FLOW_VARIABLES_H
#define GRINS_PRIMITIVE_FLOW_VARIABLES_H

// libMesh forward declarations
class GetPot;
namespace libMesh
{
  class FEMSystem;
}

// GRINS
#include "grins/variables_base.h"

namespace GRINS
{

  class PrimitiveFlowVariables : public VariablesBase
  {
  public:

    PrimitiveFlowVariables( const GetPot& input );
    ~PrimitiveFlowVariables(){};

    virtual void init( libMesh::FEMSystem* system );

    VariableIndex u_var() const;
    VariableIndex v_var() const;
    VariableIndex w_var() const;
    VariableIndex p_var() const;

  protected:

    unsigned int _u_idx, _v_idx, _w_idx, _p_idx;

  private:

    PrimitiveFlowVariables();

  };

  inline
  VariableIndex PrimitiveFlowVariables::u_var() const
  {
    return this->_vars[_u_idx];
  }

  inline
  VariableIndex PrimitiveFlowVariables::v_var() const
  {
    return this->_vars[_v_idx];
  }

  inline
  VariableIndex PrimitiveFlowVariables::w_var() const
  {
    return this->_vars[_w_idx];
  }

  inline
  VariableIndex PrimitiveFlowVariables::p_var() const
  {
    return this->_vars[_p_idx];
  }

} // end namespace GRINS

#endif //GRINS_PRIMITIVE_FLOW_VARIABLES_H
