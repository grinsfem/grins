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

#include "grins/heat_transfer_source.h"

namespace GRINS
{

  template< class SourceFunction >
  HeatTransferSource<SourceFunction>::HeatTransferSource( const std::string& physics_name, const GetPot& input )
    : HeatTransferBase(physics_name,input),
      _source(input)
  {
    return;
  }

  template< class SourceFunction >
  HeatTransferSource<SourceFunction>::~HeatTransferSource()
  {
    return;
  }

  template< class SourceFunction >
  void HeatTransferSource<SourceFunction>::element_time_derivative( bool /*compute_jacobian*/,
								    libMesh::FEMContext& context )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("HeatTransferSource::element_time_derivative");
#endif
  
    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.dof_indices_var[_T_var].size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.element_fe_var[_T_var]->get_JxW();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.element_fe_var[_T_var]->get_phi();

    // Locations of quadrature points
    const std::vector<libMesh::Point>& x_qp = context.element_fe_var[_T_var]->get_xyz();

    // Get residuals
    libMesh::DenseSubVector<Number> &FT = *context.elem_subresiduals[_T_var]; // R_{T}

    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.element_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	Real q = _source( x_qp[qp] );

	for (unsigned int i=0; i != n_T_dofs; i++)
	  {
	    FT(i) += q*T_phi[i][qp]*JxW[qp];
	  }
      }

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("HeatTransferSource::element_time_derivative");
#endif

    return;
  }

} // namespace GRINS

// Instantiate
template class GRINS::HeatTransferSource<GRINS::ConstantSourceFunction>;
