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

#include "heat_transfer_source.h"

template< class SourceFunction >
GRINS::HeatTransferSource<SourceFunction>::HeatTransferSource( const std::string& physics_name, const GetPot& input )
  : HeatTransferBase(physics_name,input),
    _source(input)
{
  this->read_input_options(input);
  return;
}

template< class SourceFunction >
GRINS::HeatTransferSource<SourceFunction>::~HeatTransferSource()
{
  return;
}

template< class SourceFunction >
void GRINS::HeatTransferSource<SourceFunction>::read_input_options( const GetPot& input )
{
  return;
}

template< class SourceFunction >
bool GRINS::HeatTransferSource<SourceFunction>::element_time_derivative( bool request_jacobian,
									 libMesh::DiffContext& context,
									 libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("HeatTransferSource::element_time_derivative");
#endif
  
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // The number of local degrees of freedom in each variable.
  const unsigned int n_T_dofs = c.dof_indices_var[_T_var].size();

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[_T_var]->get_JxW();

  // The temperature shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& T_phi =
    c.element_fe_var[_T_var]->get_phi();

  // Locations of quadrature points
  const std::vector<libMesh::Point>& x_qp = c.element_fe_var[_T_var]->get_xyz();

  // Get residuals
  libMesh::DenseSubVector<Number> &FT = *c.elem_subresiduals[_T_var]; // R_{T}

  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      Real q = _source( x_qp[qp] );

      for (unsigned int i=0; i != n_T_dofs; i++)
        {
	  FT(i) += q*T_phi[i][qp]*JxW[qp];
	}
    }

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("HeatTransferSource::element_time_derivative");
#endif

  return request_jacobian;
}

template< class SourceFunction >
bool GRINS::HeatTransferSource<SourceFunction>::side_time_derivative( bool request_jacobian,
								      libMesh::DiffContext&,
								      libMesh::FEMSystem* )
{
  return request_jacobian;
}

template< class SourceFunction >
bool GRINS::HeatTransferSource<SourceFunction>::element_constraint( bool request_jacobian,
								    libMesh::DiffContext&,
								    libMesh::FEMSystem* )
{
  return request_jacobian;
}

template< class SourceFunction >
bool GRINS::HeatTransferSource<SourceFunction>::side_constraint( bool request_jacobian,
								 libMesh::DiffContext&,
								 libMesh::FEMSystem* )
{
  return request_jacobian;
}

template< class SourceFunction >
bool GRINS::HeatTransferSource<SourceFunction>::mass_residual( bool request_jacobian,
							       libMesh::DiffContext&,
							       libMesh::FEMSystem* )
{
  return request_jacobian;
}

// Instantiate
template class GRINS::HeatTransferSource<GRINS::ConstantSourceFunction>;
