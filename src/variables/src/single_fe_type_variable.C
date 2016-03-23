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
#include "grins/single_fe_type_variable.h"

// GRINS
#include "grins/common.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"

namespace GRINS
{

  SingleFETypeVariable::SingleFETypeVariable( const GetPot& input,
                                              const std::string& physics_name,
                                              const std::string& old_var_suffix,
                                              const std::string& subsection,
                                              const std::string& default_family,
                                              const std::string& default_order,
                                              bool _is_constraint_var )
    :  FEVariablesBase(_is_constraint_var)
  {
    _family.resize(1,libMesh::INVALID_FE);
    _order.resize(1,libMesh::INVALID_ORDER);

    this->parse_family_and_order(input,
                                 physics_name,
                                 old_var_suffix,
                                 subsection,
                                 _family,
                                 _order,
                                 default_family,
                                 default_order);

    libmesh_assert_not_equal_to( _family[0], libMesh::INVALID_FE);
    libmesh_assert_not_equal_to( _order[0], libMesh::INVALID_ORDER);
  }

  SingleFETypeVariable::SingleFETypeVariable( const GetPot& input,
                                              const std::string& subsection,
                                              bool _is_constraint_var)
    :  FEVariablesBase(_is_constraint_var)
  {
     _family.resize(1,libMesh::INVALID_FE);
     _order.resize(1,libMesh::INVALID_ORDER);

     this->parse_new_style(input, subsection, _family[0], _order[0]);

     libmesh_assert_not_equal_to( _family[0], libMesh::INVALID_FE);
     libmesh_assert_not_equal_to( _order[0], libMesh::INVALID_ORDER);
  }

  void SingleFETypeVariable::parse_family_and_order( const GetPot& input,
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
            this->parse_old_style_with_warning(input,physics_name,"",default_family,default_order,
                                               subsection,family[0],order[0]);
          }
        // Some of the old style didn't have the "prefix"
        else if( input.have_variable("Physics/"+physics_name+"/FE_family") ||
                 input.have_variable("Physics/"+physics_name+"/order"))
          {
            this->parse_old_style_with_warning(input,physics_name,"",default_family,default_order,
                                               subsection,family[0],order[0]);
          }
        // We must be using the new style
        else
          {
            this->parse_new_style(input,subsection,family[0],order[0]);
          }
      }
  }

  void SingleFETypeVariable::parse_old_style_with_warning( const GetPot& input,
                                                           const std::string& physics_name,
                                                           const std::string& old_var_suffix,
                                                           const std::string& default_family,
                                                           const std::string& default_order,
                                                           const std::string& subsection,
                                                           GRINSEnums::FEFamily& family,
                                                           GRINSEnums::Order& order )
  {
    std::string warning = "WARNING: Specifying Physics/"+physics_name+"/"+old_var_suffix+"FE_family and\n";
    warning += "         Physics/"+physics_name+"/"+old_var_suffix+"order is DEPRECATED! Please update and use\n";
    warning += "         [Variables/"+subsection+"/fe_family] and [Variables/"+subsection+"/order].\n";
    grins_warning(warning);

    std::string family_in = input("Variables/"+subsection+"/fe_family", default_family );
    std::string order_in = input("Variables/"+subsection+"/order", default_order );

    family = libMesh::Utility::string_to_enum<GRINSEnums::FEFamily>(family_in);
    order = libMesh::Utility::string_to_enum<GRINSEnums::Order>(order_in);
  }

  void SingleFETypeVariable::dup_family_order_check( const GetPot& input,
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

  void SingleFETypeVariable::parse_new_style( const GetPot& input,
                                              const std::string& subsection,
                                              GRINSEnums::FEFamily& family,
                                              GRINSEnums::Order& order )
  {
    // Both options must've been specified
    if( !input.have_variable("Variables/"+subsection+"/fe_family") )
      libmesh_error_msg("ERROR: Could not find Variables/"+subsection+"/fe_family in input!");

    if( !input.have_variable("Variables/"+subsection+"/order") )
      libmesh_error_msg("ERROR: Could not find Variables/"+subsection+"/order in input!");

    std::string family_in = input("Variables/"+subsection+"/fe_family", "DIE!");
    std::string order_in = input("Variables/"+subsection+"/order", "DIE!");

    family = libMesh::Utility::string_to_enum<GRINSEnums::FEFamily>(family_in);
    order = libMesh::Utility::string_to_enum<GRINSEnums::Order>(order_in);
  }

  bool SingleFETypeVariable::have_family_or_order( const GetPot& input,
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
