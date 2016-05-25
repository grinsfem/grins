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


#ifndef GRINS_FE_VARIABLES_BASE_H
#define GRINS_FE_VARIABLES_BASE_H

// C++
#include <vector>
#include <string>
#include <limits>

// GRINS
#include "grins/grins_enums.h"
#include "grins/var_typedefs.h"
#include "grins/variables_parsing.h"

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"

// libMesh forward declarations
class GetPot;
namespace libMesh
{
  class FEMSystem;
}

namespace GRINS
{
  class FEVariablesBase
  {
  public:

    FEVariablesBase( bool is_constraint_var )
      : _is_constraint_var(is_constraint_var),
        _neumann_bc_sign(1.0)
    {}

    FEVariablesBase( const std::vector<std::string>& var_names,
                     const std::vector<VariableIndex>& var_indices )
      : _vars(var_indices),
        _var_names(var_names),
        _is_constraint_var(false),
        _neumann_bc_sign(1.0)
    {
      libmesh_assert_equal_to(var_names.size(), var_indices.size());
    }

    ~FEVariablesBase(){};

    //! Add variables to the system
    /*! This expects that _var_names has been setup during construction
        time. Most subclasses should be able to use default_fe_init, once
        they subclass this and VariablesBase. */
    virtual void init( libMesh::FEMSystem* /*system*/ ){};

    //! Set whether or not this is a "constraint" variable
    /*! Constraint variables are things like Lagrange multipliers.
        The primary implication is that if the Variable is a constraint
        variable, then no boundary conditions are required/used for
        that variable. */
    void set_is_constraint_var( bool is_constraint_var )
    { _is_constraint_var = is_constraint_var; }

    bool is_constraint_var() const
    { return _is_constraint_var; }

    //! Reset whetever Neumann bc is postive or not
    /*! Postive means a value of 1.0 will be used in front of
        NeumannBC terms while is_positive = false indicates a
        value of -1.0 should be used.*/
    void set_neumann_bc_is_positive( bool is_positive );

    libMesh::Real neumann_bc_sign() const
    { return _neumann_bc_sign; }

    //! Return the var names that are active from this class
    /*! This must not be called until init_vars has been called. */
    const std::vector<std::string>& active_var_names() const
    { return _var_names; }

    const std::vector<VariableIndex>& var_indices() const
    { return _vars; }

  protected:

    //! Method to parse variable names from input
    /*! Names parsed from: [Variables/<subsection>/names] and then
        populated into the supplied var_names vector. It is assumed
        that var_names has been properly sized, that default_names
        and var_names have the same size, and that default_names has
        been populated with unique strings. */
    void parse_names_from_input( const GetPot& input,
                                 const std::string& subsection,
                                 std::vector<std::string>& var_names,
                                 const std::vector<std::string>& default_names );

    //! Check for old name style and new name style. If both present, error.
    /*! Old name style: [Physics/VariableNames]
        New name style: [Variables/<variable type>]
        Here, we just check for the presence of the sections [Physics/VariableNames]
        and [Variables]. */
    void duplicate_name_section_check( const GetPot& input ) const;

    //! Check for deprecated variable name input style
    /*! If found, this returns true and emits a deprecated warning.
        Otherwise, this returns false.
        The string argument is supplied by each variable
        class for the warning message. E.g. if the variable class
        is going to look in "Displacement", i.e.
        [Variables/Displacement/names], then "Displacement" should be
        passed. */
    bool check_dep_name_input( const GetPot& input,
                               const std::string& new_subsection ) const;

    std::vector<VariableIndex> _vars;

    std::vector<std::string> _var_names;

    std::vector<GRINSEnums::FEFamily> _family;

    std::vector<GRINSEnums::Order> _order;

    //! Tracks whether this is a constraint variable
    /*! By constraint variable, we mean a variable that is
        effectively a Lagrange multiplier. The intended use
        case is to determine whether this variable requires
        boundary conditions to be specified (constraint variables
        do not). */
    bool _is_constraint_var;

    //! Track the sign of the Neumann BC term. Defaults to 1.0.
    /*! Depending on the Physics/Variable combination, the sign in
        front of the Neumann boundary term can change. */
    libMesh::Real _neumann_bc_sign;

  };

  inline
  void FEVariablesBase::set_neumann_bc_is_positive( bool is_positive )
  {
    if(is_positive)
      _neumann_bc_sign = 1.0;
    else
      _neumann_bc_sign = -1.0;
  }

} // end namespace GRINS

#endif // GRINS_FE_VARIABLES_BASE_H
