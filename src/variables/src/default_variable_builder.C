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
#include "grins/default_variable_builder.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/variable_factory.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  void DefaultVariableBuilder::build_variables_impl( const GetPot& input,
                                                     MultiphysicsSystem& system )
  {
    std::vector<std::string> var_sections;
    this->parse_var_sections_vector( input, var_sections );

    for( std::vector<std::string>::const_iterator var_sect = var_sections.begin();
         var_sect != var_sections.end(); ++var_sect )
      {
        // Convenience
        std::string var_section( *var_sect );

        // We construct based on the Variable type so we can name the Variables whatever we want
        // If var_type is invalid, that will be detected in the build_fe_var call, since
        // we feed the var_type to VariableFactoryAbstract::build.
        std::string var_type = this->parse_var_type( input, var_section );

        // Parse names
        std::vector<std::string> var_names;
        this->parse_var_names( input, var_type, var_section, var_names );

        // Parse FE family
        std::string fe_family = this->parse_fe_family( input, var_section );

        // Parse FE order
        std::string order = this->parse_fe_order( input, var_section );

        // Add variables to system
        std::vector<VariableIndex> var_indices;
        this->add_vars_to_system( system, var_names, fe_family, order, var_indices );

        // Build FEVariablesBase object
        SharedPtr<FEVariablesBase> fe_var = this->build_fe_var( var_type, var_names, var_indices );

        // Add to VariableWarehouse
        this->add_variable_to_warehouse( fe_var, var_section );
      }
  }

  void DefaultVariableBuilder::parse_var_names( const GetPot& input,
                                                const std::string& var_type,
                                                const std::string& var_section,
                                                std::vector<std::string>& var_names ) const
  {
    // Just in case
    var_names.clear();

    VariableFactoryAbstract::set_getpot(input);
    VariableFactoryAbstract::set_var_section(VariablesParsing::variables_section()+"/"+var_section);

    var_names = VariableFactoryAbstract::build_var_names(var_type);
  }

  std::string DefaultVariableBuilder::parse_var_option( const GetPot& input,
                                                        const std::string& var_section,
                                                        const std::string& option,
                                                        const std::string& default_val ) const
  {
    std::string input_sec = VariablesParsing::variables_section()+"/"+var_section+"/"+option;
    if(!input.have_variable(input_sec))
      libmesh_error_msg("ERROR: Could not find Variable input option "+input_sec);

    return input(input_sec,default_val);
  }

} // end namespace GRINS
