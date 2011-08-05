//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - a low Mach number Navier-Stokes Finite-Element Solver
//
// Copyright (C) 2010,2011 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include "bc_types.h"
#include "var_typedefs.h"
#include "point_func_base.h"

// libMesh stuff
#include "libmesh.h"
#include "boundary_info.h"
#include "fe_base.h"
#include "fe_interface.h"
#include "mesh.h"
#include "quadrature.h"
#include "parameters.h"
#include "string_to_enum.h"
#include "fem_context.h"
#include "fem_system.h"

namespace GRINS
{
  //! Class to hold typical boundary condition methods
  /*!
    This class holds functions to apply generic versions of
    Dirichlet and Neumann boundary conditions.
  */
  class BoundaryConditions
  {
  public:

    BoundaryConditions();
    ~BoundaryConditions();
    
    /*! Applies Dirichlet boundary conditions for a single variable using the
      penalty method. */
    void apply_dirichlet( libMesh::DiffContext &context, const bool request_jacobian,
			  const GRINS::VariableIndex var, const double value, 
			  const double penalty = 1.0e16 );

    /*! Applies Dirichlet boundary conditions for a vector of variables using the
      penalty method. */
    void apply_dirichlet( libMesh::DiffContext &context, const bool request_jacobian,
			  const std::vector<GRINS::VariableIndex>& vars,
			  const std::vector<bool>& set_vars,
			  GRINS::BasePointFuncObj* func,
			  const double penalty = 1.0e16 );

    void apply_neumann( libMesh::DiffContext &context,
			GRINS::VariableIndex var1, 
			GRINS::VariableIndex var2, double value, 
			double jacobian_value );

    /*! The idea here is to pin a variable to a particular value if there is
      a null space - e.g. pressure for IncompressibleNavierStokes. */
    void pin_value( libMesh::DiffContext &context, const bool request_jacobian,
		    const GRINS::VariableIndex var, const double value,
		    const libMesh::Point& pin_location, const double penalty = 1.0 );

    GRINS::BC_TYPES string_to_enum( const std::string bc_type );

  };
}
#endif //BOUNDARY_CONDITIONS_H
