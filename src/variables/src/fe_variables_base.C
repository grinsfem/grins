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
#include "grins/fe_variables_base.h"

// GRINS
#include "grins/common.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  void FEVariablesBase::default_fe_init( libMesh::FEMSystem* system,
                                         const std::vector<std::string>& var_names,
                                         std::vector<VariableIndex>& vars ) const
  {
    const unsigned int n_vars = var_names.size();
    vars.resize(n_vars);

    for( unsigned int v = 0; v < n_vars; v++ )
      {
        vars[v] = system->add_variable( var_names[v], this->_order[0], _family[0]);
      }
  }

  void FEVariablesBase::parse_family_and_order( const GetPot& input,
                                                const std::string& physics_name,
                                                const std::string& old_var_suffix,
                                                const std::string& subsection,
                                                std::vector<GRINSEnums::FEFamily>& family,
                                                std::vector<GRINSEnums::Order>& order,
                                                const std::string& default_family,
                                                const std::string& default_order )
  {
    libmesh_assert_equal_to( family.size(), 1 );
    libmesh_assert_equal_to( family.size(), order.size() );

    // Check if FEfamily/order set in both old style and new style.
    // Errors out if both are present.
    this->dup_family_order_check(input,physics_name,old_var_suffix);

    // If there's nothing about family or order, then we'll use the defaults.
    // This is deprecated.
    if( !this->have_family_or_order(input,physics_name,old_var_suffix,subsection) )
      {
        std::string warning = "WARNING: Could not find input for FE family and order for Variable "+subsection+".\n";
        warning += "         Using default values: FEFamily = "+default_family+", Order = "+default_order+"\n";
        warning += "         THIS IS DEPRECATED. In the future, you explicitly specify these in your input file\n";
        warning += "         using: [Variables/"+subsection+"/fe_family] and [Variables/"+subsection+"/order].\n";
        grins_warning(warning);

        family[0] = libMesh::Utility::string_to_enum<GRINSEnums::FEFamily>(default_family);
        order[0] = libMesh::Utility::string_to_enum<GRINSEnums::Order>(default_order);
      }
    // So here we know we have one or the other
    else
      {
        // If we got here, there was one or the other, but it's possible they only specified one of them
        // so we default to the default value in this case.
        // Using this style of input is deprecated.
        if( input.have_variable("Physics/"+physics_name+"/"+old_var_suffix+"FE_family") ||
            input.have_variable("Physics/"+physics_name+"/"+old_var_suffix+"order"))
          {
            std::string family_in, order_in;
            this->parse_old_style_with_warning(input,physics_name,old_var_suffix,default_family,default_order,
                                               subsection,family_in,order_in);

            family[0] = libMesh::Utility::string_to_enum<GRINSEnums::FEFamily>(family_in);
            order[0] = libMesh::Utility::string_to_enum<GRINSEnums::Order>(order_in);
          }
        // Some of the old style didn't have the "prefix"
        else if( input.have_variable("Physics/"+physics_name+"/FE_family") ||
                 input.have_variable("Physics/"+physics_name+"/order"))
          {
            std::string family_in, order_in;
            this->parse_old_style_with_warning(input,physics_name,"",default_family,default_order,
                                               subsection,family_in,order_in);

            family[0] = libMesh::Utility::string_to_enum<GRINSEnums::FEFamily>(family_in);
            order[0] = libMesh::Utility::string_to_enum<GRINSEnums::Order>(order_in);
          }
        // We must be using the new style
        else
          {
            // Both options must've been specified
            if( !input.have_variable("Variables/"+subsection+"/fe_family") )
              libmesh_error_msg("ERROR: Could not find Variables/"+subsection+"/fe_family in input!");

            if( !input.have_variable("Variables/"+subsection+"/order") )
              libmesh_error_msg("ERROR: Could not find Variables/"+subsection+"/order in input!");

            std::string family_in = input("Variables/"+subsection+"/fe_family", "DIE!");
            std::string order_in = input("Variables/"+subsection+"/order", "DIE!");

            family[0] = libMesh::Utility::string_to_enum<GRINSEnums::FEFamily>(family_in);
            order[0] = libMesh::Utility::string_to_enum<GRINSEnums::Order>(order_in);
          }
      }
  }

  void FEVariablesBase::parse_old_style_with_warning( const GetPot& input,
                                                      const std::string& physics_name,
                                                      const std::string& old_var_suffix,
                                                      const std::string& default_family,
                                                      const std::string& default_order,
                                                      const std::string& subsection,
                                                      std::string& parsed_family,
                                                      std::string& parsed_order )
  {
     std::string warning = "WARNING: Specifying Physics/"+physics_name+"/"+old_var_suffix+"FE_family and\n";
     warning += "         Physics/"+physics_name+"/"+old_var_suffix+"order is DEPRECATED! Please update and use\n";
     warning += "         [Variables/"+subsection+"/fe_family] and [Variables/"+subsection+"/order].\n";
     grins_warning(warning);

     parsed_family = input("Physics/"+physics_name+"/"+old_var_suffix+"FE_family", default_family);
     parsed_order = input("Physics/"+physics_name+"/"+old_var_suffix+"order", default_order);
  }

  void FEVariablesBase::dup_family_order_check( const GetPot& input,
                                                const std::string& physics_name,
                                                const std::string& old_var_suffix ) const
  {
    // We'll error out if there's even a [Variables] section and the old style
    // in order to be conservative.
    if( ( input.have_variable("Physics/"+physics_name+"/FE_family") ||
          input.have_variable("Physics/"+physics_name+"/order") ||
          input.have_variable("Physics/"+physics_name+"/"+old_var_suffix+"FE_family") ||
          input.have_variable("Physics/"+physics_name+"/"+old_var_suffix+"order")) &&
        input.have_section("Variables") )
      libmesh_error_msg("ERROR: Cannot have a [Variables] section and deprecated FE_family/order input style!");
  }

  bool FEVariablesBase::have_family_or_order( const GetPot& input,
                                              const std::string& physics_name,
                                              const std::string& old_var_suffix,
                                              const std::string& subsection ) const
  {
    bool have_family_or_order = false;

    if( input.have_variable("Physics/"+physics_name+"/FE_family") ||
        input.have_variable("Physics/"+physics_name+"/order") ||
        input.have_variable("Physics/"+physics_name+"/"+old_var_suffix+"FE_family") ||
        input.have_variable("Physics/"+physics_name+"/"+old_var_suffix+"order") ||
        input.have_variable("Variables/"+subsection+"/fe_family") ||
        input.have_variable("Variables/"+subsection+"/order") )
      have_family_or_order = true;

    return have_family_or_order;
  }

} // end namespace GRINS
