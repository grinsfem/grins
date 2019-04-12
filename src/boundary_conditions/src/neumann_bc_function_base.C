//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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

// This class
#include "grins/neumann_bc_function_base.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/parsed_function_traits.h"

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/point.h"
#include "libmesh/quadrature.h"
#include "libmesh/elem.h"

namespace GRINS
{
  template<typename FunctionType, typename FEShape>
  bool NeumannBCFunctionBase<FunctionType,FEShape>::eval_flux( bool compute_jacobian,
                                                               AssemblyContext& context,
                                                               libMesh::Real sign,
                                                               bool is_axisymmetric )
  {
    libmesh_assert(!_vars.empty());

    unsigned int n_vars = this->_vars.size();

    // The number of local degrees of freedom in each variable.
    // We're assuming that each var in our set is the same FEType
    const unsigned int n_var_dofs = context.get_dof_indices(this->_vars[0]).size();
#ifndef NDEBUG
    for( unsigned int v = 1; v < n_vars; v++ )
      libmesh_assert_equal_to( n_var_dofs, context.get_dof_indices(this->_vars[v]).size() );
#endif

    libMesh::FEGenericBase<FEShape>* side_fe = NULL;
    context.get_side_fe( this->_vars[0], side_fe );

    // Element Jacobian * quadrature weight for side integration.
    const std::vector<libMesh::Real> &JxW_side = side_fe->get_JxW();

    // The var shape functions at side quadrature points.
    const std::vector<std::vector<FEShape> >& var_phi_side =
      side_fe->get_phi();

    // Physical location of the quadrature points
    const std::vector<libMesh::Point>& var_qpoint =
      side_fe->get_xyz();

    std::vector<libMesh::DenseSubVector<libMesh::Number>*> F_vars(_vars.size(),NULL);

    for( unsigned int v = 0; v < n_vars; v++ )
      F_vars[v] = &(context.get_elem_residual(this->_vars[v])); // residual

    unsigned int n_qpoints = context.get_side_qrule().n_points();

    FEShape value = 0.0;

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Real jac = JxW_side[qp];

        if( is_axisymmetric )
          {
            const libMesh::Number r = var_qpoint[qp](0);
            jac *= r;
          }

        for( unsigned int v = 0; v < n_vars; v++ )
          {
            value = this->eval_func( context, var_qpoint[qp], context.get_time(), this->_vars[v], *_func );

            for (unsigned int i=0; i != n_var_dofs; i++)
              (*F_vars[v])(i) += sign*(value*var_phi_side[i][qp])*jac;
          }
      }

    if( ParsedFunctionTraits<FunctionType>::is_fem_function && compute_jacobian )
      libmesh_error_msg("ERROR: Analytical Jacobians of FEMFunctionBase objects currently not supported!");

    return compute_jacobian;
  }

  // Instantiate
  template class NeumannBCFunctionBase<libMesh::FunctionBase<libMesh::Number> ,libMesh::Number>;
  template class NeumannBCFunctionBase<libMesh::FEMFunctionBase<libMesh::Number> ,libMesh::Number>;
  template class NeumannBCFunctionBase<libMesh::FunctionBase<libMesh::Gradient> ,libMesh::Gradient>;
  template class NeumannBCFunctionBase<libMesh::FEMFunctionBase<libMesh::Gradient> ,libMesh::Gradient>;

} // end namespace GRINS
