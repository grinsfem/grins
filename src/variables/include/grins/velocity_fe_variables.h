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


#ifndef GRINS_VELOCITY_FE_VARIABLES_H
#define GRINS_VELOCITY_FE_VARIABLES_H

// GRINS
#include "grins/multi_var_single_fe_type_variable.h"
#include "grins/variables_parsing.h"

namespace GRINS
{

  class VelocityFEVariables : public MultiVarSingleFETypeVariable
  {
  public:

    VelocityFEVariables( const GetPot& input, const std::string& physics_name,
                         bool is_constraint_var = false );
    ~VelocityFEVariables(){};

    virtual void init( libMesh::FEMSystem* system );

    VariableIndex u() const;
    VariableIndex v() const;
    VariableIndex w() const;

  private:

    VelocityFEVariables();

    std::string subsection() const
    { return VariablesParsing::velocity_section(); }

    std::vector<std::string> old_var_names()
    {
      std::vector<std::string> var_names(3);
      var_names[0] = "u_velocity";
      var_names[1] = "v_velocity";
      var_names[2] = "w_velocity";
      return var_names;
    }

    std::vector<std::string> default_names()
    {
      std::vector<std::string> var_names(3);
      var_names[0] = "u";
      var_names[1] = "v";
      var_names[2] = "w";
      return var_names;
    }

    unsigned int _u_idx, _v_idx, _w_idx;

  };

  inline
  VariableIndex VelocityFEVariables::u() const
  {
    return this->_vars[_u_idx];
  }

  inline
  VariableIndex VelocityFEVariables::v() const
  {
    return this->_vars[_v_idx];
  }

  inline
  VariableIndex VelocityFEVariables::w() const
  {
    return this->_vars[_w_idx];
  }

} // end namespace GRINS

#endif //GRINS_VELOCITY_FE_VARIABLES_H
