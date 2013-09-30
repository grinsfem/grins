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

#ifndef GRINS_BOUNDARY_CONDITIONS_H
#define GRINS_BOUNDARY_CONDITIONS_H

// GRINS
#include "grins/var_typedefs.h"
#include "grins/neumann_func_obj.h"

// libMesh
#include "libmesh/libmesh.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_tools.h"

namespace GRINS
{
  // Forward Declarations
  class AssemblyContext;
  class CachedValues;

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
    template<typename FEShape = libMesh::Real>
    void apply_neumann( AssemblyContext& context,
			const VariableIndex var,
			const libMesh::Real sign,
			const typename libMesh::TensorTools::IncrementRank<FEShape>::type& value ) const;

    //! Applies Neumann boundary conditions for the constant case.
    template<typename FEShape = libMesh::Real>
    void apply_neumann_axisymmetric( AssemblyContext& context,
				     const VariableIndex var,
				     const libMesh::Real sign,
				     const typename libMesh::TensorTools::IncrementRank<FEShape>::type& value ) const;

    //! Applies Neumann boundary conditions for the constant case.
    /*! This method is for the case where Neumann boundary condition is
        not in terms of a flux vector, but rather only the normal component.*/
    template<typename FEShape = libMesh::Real>
    void apply_neumann_normal( AssemblyContext& context,
			       const VariableIndex var,
			       const libMesh::Real sign,
			       const FEShape& value ) const;

    //! Applies Neumann boundary conditions for the constant case.
    /*! This method is for the case where Neumann boundary condition is
        not in terms of a flux vector, but rather only the normal component.*/
    template<typename FEShape = libMesh::Real>
    void apply_neumann_normal_axisymmetric( AssemblyContext& context,
                                            const VariableIndex var,
                                            const libMesh::Real sign,
                                            const FEShape& value ) const;


    //! Applies Neumann boundary conditions using a user-supplied function.
    /*! This function must also be aware of the Jacobian with respect to other variables. */
    template<typename FEShape = libMesh::Real>
    void apply_neumann( AssemblyContext& context,
			const CachedValues& cache,
			const bool request_jacobian,
			const VariableIndex var,
			const libMesh::Real sign,
			std::tr1::shared_ptr<NeumannFuncObj> neumann_func  ) const;

    //! Applies Neumann boundary conditions using a user-supplied function.
    /*! This function must also be aware of the Jacobian with respect to other variables. */
    template<typename FEShape = libMesh::Real>
    void apply_neumann_axisymmetric( AssemblyContext& context,
				     const CachedValues& cache,
				     const bool request_jacobian,
				     const VariableIndex var,
				     const libMesh::Real sign,
				     const std::tr1::shared_ptr<NeumannFuncObj> neumann_func  ) const;

    //! Applies Neumann boundary conditions using a user-supplied function.
    /*!  This method is for the case where Neumann boundary condition is
         not in terms of a flux vector, but rather only the normal component.
         This function must also be aware of the Jacobian with respect to other variables. */
    template<typename FEShape = libMesh::Real>
    void apply_neumann_normal( AssemblyContext& context,
			       const CachedValues& cache,
			       const bool request_jacobian,
			       const VariableIndex var,
			       const libMesh::Real sign,
			       std::tr1::shared_ptr<NeumannFuncObj> neumann_func  ) const;

    //! Applies Neumann boundary conditions using a user-supplied function.
    /*!  This method is for the case where Neumann boundary condition is
         not in terms of a flux vector, but rather only the normal component.
         This function must also be aware of the Jacobian with respect to other variables. */
    template<typename FEShape = libMesh::Real>
    void apply_neumann_normal_axisymmetric( AssemblyContext& context,
                                            const CachedValues& cache,
                                            const bool request_jacobian,
                                            const VariableIndex var,
                                            const libMesh::Real sign,
                                            std::tr1::shared_ptr<NeumannFuncObj> neumann_func  ) const;

    /*! The idea here is to pin a variable to a particular value if there is
      a null space - e.g. pressure for IncompressibleNavierStokes. */
    template<typename FEShape = libMesh::Real>
    void pin_value( AssemblyContext& context, const CachedValues& cache,
		    const bool request_jacobian,
		    const VariableIndex var, const double value,
		    const libMesh::Point& pin_location, const double penalty = 1.0 );

  };

} // end namespace GRINS
#endif // GRINS_BOUNDARY_CONDITIONS_H
