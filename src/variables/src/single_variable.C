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

// This class
#include "grins/single_variable.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  SingleVariable::SingleVariable( const GetPot& input,
                                  const std::string& old_var_name,
                                  const std::string& subsection,
                                  const std::string& default_name )
  {
    _vars.resize(1,invalid_var_index);
    _var_names.resize(1);

    std::vector<std::string> default_names(1,default_name);

    if( this->check_dep_name_input(input,subsection) )
      _var_names[0] = input("Physics/VariableNames/"+old_var_name, default_names[0] );
    else
      this->parse_names_from_input(input,subsection,_var_names,default_names);
  }

  SingleVariable::SingleVariable( const GetPot& input,
                                  const std::string& subsection,
                                  const std::string& default_name )
  {
    _vars.resize(1,invalid_var_index);
    _var_names.resize(1);

    std::vector<std::string> default_names(1,default_name);

    this->parse_names_from_input(input,subsection,_var_names,default_names);
  }

} // end namespace GRINS
