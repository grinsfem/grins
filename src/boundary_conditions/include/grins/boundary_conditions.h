//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

// GRINS
#include "grins/var_typedefs.h"
#include "grins/neumann_func_obj.h"

// libMesh
#include "libmesh/libmesh.h"

// libMesh forward declarations
namespace libMesh
{
  class Point;
  class DiffContext;
}

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

    //! Applies Neumann boundary conditions for the constant case.
    void apply_neumann( libMesh::DiffContext &context,
			const GRINS::VariableIndex var,
			const libMesh::Real sign,
			const libMesh::Point& value ) const;

    //! Applies Neumann boundary conditions using a user-supplied function.
    /*! This function must also be aware of the Jacobian with respect to other variables. */
    void apply_neumann_axisymmetric( libMesh::DiffContext &context,
				     const bool request_jacobian,
				     const GRINS::VariableIndex var,
				     const libMesh::Real sign,
				     const std::tr1::shared_ptr<GRINS::NeumannFuncObj> neumann_func  ) const;

    //! Applies Neumann boundary conditions for the constant case.
    void apply_neumann_axisymmetric( libMesh::DiffContext &context,
				     const GRINS::VariableIndex var,
				     const libMesh::Real sign,
				     const libMesh::Point& value ) const;

    //! Applies Neumann boundary conditions using a user-supplied function.
    /*! This function must also be aware of the Jacobian with respect to other variables. */
    void apply_neumann( libMesh::DiffContext &context,
			const bool request_jacobian,
			const GRINS::VariableIndex var,
			const libMesh::Real sign,
			std::tr1::shared_ptr<GRINS::NeumannFuncObj> neumann_func  ) const;

    /*! The idea here is to pin a variable to a particular value if there is
      a null space - e.g. pressure for IncompressibleNavierStokes. */
    void pin_value( libMesh::DiffContext &context, const bool request_jacobian,
		    const GRINS::VariableIndex var, const double value,
		    const libMesh::Point& pin_location, const double penalty = 1.0 );

  };
}
#endif //BOUNDARY_CONDITIONS_H
