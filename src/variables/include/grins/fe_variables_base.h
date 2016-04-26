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
#include <limits>

// GRINS
#include "grins/grins_enums.h"
#include "grins/var_typedefs.h"

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

    ~FEVariablesBase(){};

    //! Add variables to the system
    /*! This expects that _var_names has been setup during construction
        time. Most subclasses should be able to use default_fe_init, once
        they subclass this and VariablesBase. */
    virtual void init( libMesh::FEMSystem* system ) =0;

    bool is_constraint_var() const
    { return _is_constraint_var; }

    static void set_is_axisymmetric( bool is_axisymmetric )
    { _is_axisymmetric = is_axisymmetric; }

    static bool is_axisymmetric()
    { return _is_axisymmetric; }

    //! Reset Neumann bc sign to 1.0 or -1.0.
    /*! Error is thrown if incoming neumann_bc_sign does not have
        magnitude 1. */
    void reset_neumann_bc_sign( libMesh::Real neumann_bc_sign );

    libMesh::Real neumann_bc_sign() const
    { return _neumann_bc_sign; }

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

    //! Track whether this is an axisymmetric problem
    /*! This is static because either everyone thinks the problem
        axisymmetric or not. */
    static bool _is_axisymmetric;

    //! Track the sign of the Neumann BC term. Defaults to 1.0.
    /*! Depending on the Physics/Variable combination, the sign in
        front of the Neumann boundary term can change. */
    libMesh::Real _neumann_bc_sign;

  };

  inline
  void FEVariablesBase::reset_neumann_bc_sign( libMesh::Real neumann_bc_sign )
  {
    _neumann_bc_sign = neumann_bc_sign;
    if( std::abs( _neumann_bc_sign - 1.0 ) > std::numeric_limits<libMesh::Real>::epsilon() )
      libmesh_error_msg("ERROR: neumann_bc_sign must be 1.0 or -1.0!");
  }

} // end namespace GRINS

#endif // GRINS_FE_VARIABLES_BASE_H
