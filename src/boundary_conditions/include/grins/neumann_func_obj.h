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
#ifndef GRINS_NEUMANN_FUNC_OBJ_H 
#define GRINS_NEUMANN_FUNC_OBJ_H

// GRINS
#include "grins/var_typedefs.h"

// libMesh
#include "libmesh/point.h"

// C++
#include <vector>

// libMesh forward declarations
namespace libMesh
{
  class FEMContext;
}

namespace GRINS
{
  // GRINS forward declarations
  class CachedValues;

  //! Base class for general, non-constant Neumann boundary conditions
  class NeumannFuncObj
  {
  public:
    
    NeumannFuncObj();
    
    virtual ~NeumannFuncObj();
    
    //! Returns the value of the implemented Neumann boundary condition
    /*! This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp. */
    virtual libMesh::Point value( const libMesh::FEMContext& context,
				  const CachedValues& cache, const unsigned int qp );

    //! Returns the value of the implemented Neumann boundary condition
    /*! This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp. 
      Returns the normal component of the Neumann value. Only to be used when flux vector is
      formulated implicitly in terms of normal component. By default, does nothing since
      it's only applicable in special cases. */
    virtual libMesh::Real normal_value( const libMesh::FEMContext& context, const CachedValues& cache,
					const unsigned int qp );
    
    //! Returns the derivative with respect to the primary variable of the implemented Neumann boundary condition.
    /*! This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp. */
    virtual libMesh::Point derivative( const libMesh::FEMContext& context, const CachedValues& cache,
				       const unsigned int qp );

    //! Returns the derivative with respect to the primary variable of the implemented Neumann boundary condition.
    /*! This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp.
      Returns the normal component of the Neumann value. Only to be used when flux vector is
      formulated implicitly in terms of normal component. By default, does nothing since
      it's only applicable in special cases. */
    virtual libMesh::Real normal_derivative( const libMesh::FEMContext& context, const CachedValues& cache,
					     const unsigned int qp );

    //! If needed, returns the derivative with respect to other variables in the system.
    /*! By default, does nothing. User should reimplement is this is needed.
      This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp. */
    virtual libMesh::Point derivative( const libMesh::FEMContext& context, const CachedValues& cache,
				       const unsigned int qp, 
				       const GRINS::VariableIndex jac_var );

    //! If needed, returns the derivative with respect to other variables in the system.
    /*! By default, does nothing. User should reimplement is this is needed.
      This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp.
      Returns the normal component of the Neumann value. Only to be used when flux vector is
      formulated implicitly in terms of normal component. By default, does nothing since
      it's only applicable in special cases. */
    virtual libMesh::Real normal_derivative( const libMesh::FEMContext& context, const CachedValues& cache,
					     const unsigned int qp, 
					     const GRINS::VariableIndex jac_var );

    const std::vector<VariableIndex>& get_other_jac_vars();

  protected:

    std::vector<VariableIndex> _jac_vars;

  }; // class NeumannFuncObj
  
} // end namespace GRINS
#endif // GRINS_NEUMANN_FUNC_OBJ_H
