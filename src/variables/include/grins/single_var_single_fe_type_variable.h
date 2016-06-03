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

#ifndef GRINS_SINGLE_VAR_SINGLE_FE_TYPE_VARIABLE_H
#define GRINS_SINGLE_VAR_SINGLE_FE_TYPE_VARIABLE_H

// GRINS
#include "grins/single_fe_type_variable.h"

namespace GRINS
{
  //! Base class for variables with a single component
  class SingleVarSingleFETypeVariable : public SingleFETypeVariable
  {
  public:

    //! Deprecated, old style constructor
    /*! This constructor is used for when there is possibly old deprecated
        styles of input for which we do additional checks/warnings. Otherwise,
        you should use the new constructor. */
    SingleVarSingleFETypeVariable( const GetPot& input,
                                   const std::string& physics_name,
                                   const std::string& old_var_prefix,
                                   const std::string& old_var_name,
                                   const std::string& default_name,
                                   const std::string& subsection,
                                   const std::string& default_family,
                                   const std::string& default_order,
                                   bool is_constraint_var );

    //! Primary constructor
    /*! Will parse from input section [Variables/<subsection>]. */
    SingleVarSingleFETypeVariable( const GetPot& input,
                                   const std::string& subsection,
                                   const std::string& default_name,
                                   bool is_constraint_var);

    ~SingleVarSingleFETypeVariable(){};

  };

} // end namespace GRINS

#endif // GRINS_SINGLE_VAR_SINGLE_FE_TYPE_VARIABLE_H
