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

#ifndef GRINS_TURBULENCE_FE_VARIABLES_H
#define GRINS_TURBULENCE_FE_VARIABLES_H

// GRINS
#include "grins/single_var_single_fe_type_variable.h"
#include "grins/variables_parsing.h"

// libMesh forward declarations
class GetPot;
namespace libMesh
{
  class FEMSystem;
}

namespace GRINS
{
  class TurbulenceFEVariables : public SingleVarSingleFETypeVariable
  {
  public:

    TurbulenceFEVariables( const GetPot& input, const std::string& physics_name,
                           bool _is_constraint_var = false )
      :  SingleVarSingleFETypeVariable(input,physics_name,"TU_",this->old_var_name(),this->default_name(),
                                       this->subsection(),"LAGRANGE","FIRST",_is_constraint_var)
    {}

    TurbulenceFEVariables( const std::vector<std::string>& var_names,
                           const std::vector<VariableIndex>& var_indices )
      : SingleVarSingleFETypeVariable(var_names,var_indices)
    {}

    ~TurbulenceFEVariables(){};

    VariableIndex nu() const;

  private:

    TurbulenceFEVariables();

    std::string old_var_name() const
    { return "turbulent_viscosity"; }

    std::string subsection() const
    { return VariablesParsing::turbulence_section(); }

    std::string default_name() const
    { return "nu"; }

  };

  inline
  VariableIndex TurbulenceFEVariables::nu() const
  {
    return _vars[0];
  }

} // end namespace GRINS

#endif // GRINS_TURBULENCE_FE_VARIABLES_H
