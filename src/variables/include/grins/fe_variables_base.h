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


#ifndef GRINS_FE_VARIABLES_BASE_H
#define GRINS_FE_VARIABLES_BASE_H

// C++
#include <vector>
#include <string>

// GRINS
#include "grins/grins_enums.h"
#include "grins/var_typedefs.h"

// libMesh
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
      : _is_constraint_var(is_constraint_var)
    {}

    ~FEVariablesBase(){};

    //! Add variables to the system
    /*! This expects that _var_names has been setup during construction
        time. Most subclasses should be able to use default_fe_init, once
        they subclass this and VariablesBase. */
    virtual void init( libMesh::FEMSystem* system ) =0;

    bool is_constraint_var() const
    { return _is_constraint_var; }

  protected:

    //! Default method for init'ing variables
    /*! This method will init all variables in the var_names
        vector using the first _family and the first _order
        finite elements and populate the vars vector with the
        corresponding variable numbers. If the user has more than
        one FE type, they should not use this method. */
    void default_fe_init( libMesh::FEMSystem* system,
                          const std::vector<std::string>& var_names,
                          std::vector<VariableIndex>& vars ) const;

    std::vector<GRINSEnums::FEFamily> _family;

    std::vector<GRINSEnums::Order> _order;

    //! Tracks whether this is a constraint variable
    /*! By constraint variable, we mean a variable that is
        effectively a Lagrange multiplier. The intended use
        case is to determine whether this variable requires
        boundary conditions to be specified (constraint variables
        do not). This should be set by the finite element
        type subclasses. */
    bool _is_constraint_var;

  };
} // end namespace GRINS

#endif // GRINS_FE_VARIABLES_BASE_H
