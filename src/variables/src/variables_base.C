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
#include "grins/variables_base.h"

// GRINS
#include "grins/common.h"

// libMesh
#include "libmesh/fem_system.h"
#include "libmesh/getpot.h"

namespace GRINS
{
  void VariablesBase::default_var_init( libMesh::FEMSystem* system )
  {
    unsigned int n_vars = _var_names.size();

    _vars.resize(n_vars);

    for( unsigned int v = 0; v < n_vars; v++ )
      {
        libmesh_assert( system->has_variable(_var_names[v]) );
        _vars[v] = system->variable_number(_var_names[v]);
      }
  }

  void VariablesBase::parse_names_from_input( const GetPot& input,
                                              const std::string& subsection,
                                              std::vector<std::string>& var_names,
                                              const std::vector<std::string>& default_names )
  {
    libmesh_assert_equal_to( var_names.size(), default_names.size() );

    unsigned int n_names = default_names.size();

    for( unsigned int n = 0; n < n_names; n++ )
      var_names[n] = input("Variables/"+subsection+"/names", default_names[n], n);
  }

  void VariablesBase::duplicate_name_section_check( const GetPot& input ) const
  {
    if( input.have_section("Physics/VariableNames") &&
        input.have_section("Variables") )
      libmesh_error_msg("ERROR: Cannot have both Physics/VariableNames and Variables in input!");
  }

  bool VariablesBase::check_dep_name_input( const GetPot& input,
                                            const std::string& new_subsection ) const
  {
    this->duplicate_name_section_check(input);

    bool is_old_input_style = false;

    if( input.have_section("Physics/VariableNames") )
      {
        is_old_input_style = true;

        std::string warning = "WARNING: Specifying variable names with Physics/VariableNames is DEPRECATED!\n";
        warning += "         Please update to use Variables/"+new_subsection+"/names.\n";
        grins_warning(warning);
      }

    return is_old_input_style;
  }

} // end namespace GRINS
