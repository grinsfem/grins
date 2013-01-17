//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
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
#ifndef NEUMANN_FUNC_OBJ_H 
#define NEUMANN_FUNC_OBJ_H

// libMesh stuff
#include "libmesh.h"
#include "point.h"
#include "fem_context.h"

// GRINS stuff
#include "grins/var_typedefs.h"

namespace GRINS
{

  //! Abstract base class for general, non-constant Neumann boundary conditions
  class NeumannFuncObj
  {
  public:
    
    NeumannFuncObj();
    
    virtual ~NeumannFuncObj();
    
    //! Returns the value of the implemented Neumann boundary condition
    /*! This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp. */
    virtual libMesh::Point value( const libMesh::FEMContext& context, const unsigned int qp ) = 0;
    
    //! Returns the derivative with respect to the primary variable of the implemented Neumann boundary condition.
    /*! This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp. */
    virtual libMesh::Point derivative( const libMesh::FEMContext& context, const unsigned qp ) = 0;

    //! If needed, returns the derivative with respect to other variables in the system.
    /*! By default, does nothing. User should reimplement is this is needed.
      This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp. */
    virtual libMesh::Point derivative( const libMesh::FEMContext& context, const unsigned int qp, 
				       const GRINS::VariableIndex jac_var ); 

    std::vector<GRINS::VariableIndex> get_other_jac_vars();

  protected:

    std::vector<GRINS::VariableIndex> _jac_vars;

  }; // class NeumannFuncObj
  
} // namespace GRINS

#endif // NEUMANN_FUNC_OBJ_H
