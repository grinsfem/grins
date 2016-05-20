//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_THERMO_PRESSURE_FE_VARIABLE_H
#define GRINS_THERMO_PRESSURE_FE_VARIABLE_H

// GRINS
#include "grins/single_var_single_fe_type_variable.h"
#include "grins/variables_parsing.h"

namespace GRINS
{

  class ThermoPressureFEVariable : public SingleVarSingleFETypeVariable
  {
  public:

    ThermoPressureFEVariable( const GetPot& input, const std::string& physics_name,
                              bool _is_constraint_var = false )
      :  SingleVarSingleFETypeVariable(input,physics_name,"",this->old_var_name(),this->default_name(),
                                       this->subsection(),"SCALAR","FIRST",_is_constraint_var)
    {
      // Currently only support SCALAR and FIRST
      libmesh_assert_equal_to( _family[0], libMesh::SCALAR );
      libmesh_assert_equal_to( _order[0], libMesh::FIRST );
    }

    ThermoPressureFEVariable( const std::vector<std::string>& var_names,
                              const std::vector<VariableIndex>& var_indices )
      : SingleVarSingleFETypeVariable(var_names,var_indices)
    {}

    virtual ~ThermoPressureFEVariable(){};

    VariableIndex p0() const;

  private:

    ThermoPressureFEVariable();

    std::string old_var_name() const
    { return "thermo_presure"; }

    std::string subsection() const
    { return VariablesParsing::thermo_pressure_section(); }

    std::string default_name() const
    { return "p0"; }

  };

  inline
  VariableIndex ThermoPressureFEVariable::p0() const
  {
    return _vars[0];
  }

} // end namespace GRINS

#endif // GRINS_THERMO_PRESSURE_FE_VARIABLE_H
