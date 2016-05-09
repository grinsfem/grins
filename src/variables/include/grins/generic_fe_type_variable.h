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

#ifndef GRINS_GENERIC_FE_TYPE_VARIABLE_H
#define GRINS_GENERIC_FE_TYPE_VARIABLE_H

// libMesh forward declarations
class GetPot;

// GRINS
#include "grins/single_var_single_fe_type_variable.h"
#include "grins/variables_parsing.h"

namespace GRINS
{
  //! Generic FE variable for generic physics
  /*! The variable inputs, e.g. fe_family, will be tied to the input physics_name.
      Thus, the input specification will be [Variables/GenericVariable:<physics_name>/fe_family],
      etc.*/
  class GenericFETypeVariable : public SingleVarSingleFETypeVariable
  {
  public:

    GenericFETypeVariable( const GetPot& input,
                           const std::string& physics_name,
                           bool _is_constraint_var = false )
      : SingleVarSingleFETypeVariable(input,this->section_name(physics_name),this->default_name(),_is_constraint_var)
    {}

    ~GenericFETypeVariable(){};

    VariableIndex var() const;

  protected:

    GenericFETypeVariable();

    std::string section_name(const std::string& physics_name) const
    { return VariablesParsing::generic_section()+":"+physics_name; }

    std::string default_name() const
    { return "u"; }

  };

  inline
  VariableIndex GenericFETypeVariable::var() const
  {
    return _vars[0];
  }

} // end namespace GRINS

#endif // GRINS_GENERIC_FE_TYPE_VARIABLE_H
