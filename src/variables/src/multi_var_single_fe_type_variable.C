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

// This class
#include "grins/multi_var_single_fe_type_variable.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  MultiVarSingleFETypeVariable::MultiVarSingleFETypeVariable( const GetPot& input,
                                                              const std::string& physics_name,
                                                              const std::string& old_var_prefix,
                                                              const std::vector<std::string>& old_var_names,
                                                              const std::vector<std::string>& default_names,
                                                              const std::string& subsection,
                                                              const std::string& default_family,
                                                              const std::string& default_order,
                                                              bool is_constraint_var )
  : SingleFETypeVariable(input,physics_name,old_var_prefix,subsection,default_family,
                         default_order,is_constraint_var)
  {
    libmesh_assert_equal_to( default_names.size(), old_var_names.size() );

    _var_names.resize(default_names.size());

    if( this->check_dep_name_input(input,subsection) )
      for( unsigned int n = 0; n < default_names.size(); n++ )
        _var_names[n] = input("Physics/VariableNames/"+old_var_names[n], default_names[n] );
    else
      this->parse_names_from_input(input,subsection,_var_names,default_names);
  }

  MultiVarSingleFETypeVariable::MultiVarSingleFETypeVariable( const GetPot& input,
                                                              const std::string& subsection,
                                                              const std::vector<std::string>& default_names,
                                                              bool is_constraint_var )
    : SingleFETypeVariable(input,subsection,is_constraint_var)
  {
    _var_names.resize(default_names.size());
    this->parse_names_from_input(input,subsection,_var_names,default_names);
  }

} // end namespace GRINS
