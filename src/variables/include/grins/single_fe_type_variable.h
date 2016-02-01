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

#ifndef GRINS_SINGLE_FE_TYPE_VARIABLE_H
#define GRINS_SINGLE_FE_TYPE_VARIABLE_H

// GRINS
#include "grins/fe_variables_base.h"

namespace GRINS
{
  class SingleFETypeVariable : public FEVariablesBase
  {
  public:

    SingleFETypeVariable( const GetPot& input,
                          const std::string& physics_name,
                          const std::string& old_var_suffix,
                          const std::string& subsection,
                          const std::string& default_family,
                          const std::string& default_order );

    ~SingleFETypeVariable(){};

  protected:

    //! Helper function to parse FEFamily and Order.
    /*! Mainly to encapsulate warning/using old style and new style.
        Note that default_family and default_order are *only* for the
        old style. In the new style, the user *must* specify the FEFamily
        and order in the input file. This function assumes that family and order
        have been properly sized. Currently assumes only a single
        family and order. */
    void parse_family_and_order( const GetPot& input,
                                 const std::string& physics_name,
                                 const std::string& old_var_suffix,
                                 const std::string& subsection,
                                 std::vector<GRINSEnums::FEFamily>& family,
                                 std::vector<GRINSEnums::Order>& order,
                                 const std::string& default_family,
                                 const std::string& default_order );

    //! Check (and error if true) for old and new style FEFamily/Order input
    /*! Actually, for the new style, we just check for the presence
        of a [Variables] section in order to be conservative. */
    void dup_family_order_check( const GetPot& input,
                                 const std::string& physics_name,
                                 const std::string& old_var_suffix) const;

    //! Check for *no* presence of FEFamily/Order input
    bool have_family_or_order( const GetPot& input,
                               const std::string& physics_name,
                               const std::string& old_var_suffix,
                               const std::string& subsection ) const;

    void parse_old_style_with_warning( const GetPot& input,
                                       const std::string& physics_name,
                                       const std::string& old_var_suffix,
                                       const std::string& default_family,
                                       const std::string& default_order,
                                       const std::string& subsection,
                                       std::string& parsed_family,
                                       std::string& parsed_order );

  };

} // end namespace GRINS

#endif // GRINS_SINGLE_FE_TYPE_VARIABLE_H
