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

// This class
#include "grins/heat_transfer_source.h"

// GRINS
#include "grins/constant_source_func.h"

// libMesh
#include "libmesh/utility.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/fem_context.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace GRINS
{

  template< class SourceFunction >
  HeatTransferSource<SourceFunction>::HeatTransferSource( const std::string& physics_name, const GetPot& input )
    : Physics(physics_name,input),
      _source(input),
      _T_var_name( input("Physics/VariableNames/Temperature", T_var_name_default ) ),
      _T_FE_family( libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+heat_transfer+"/FE_family", "LAGRANGE") ) ),
      _T_order( libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+heat_transfer+"/T_order", "SECOND") ) )
  {
    return;
  }

  template< class SourceFunction >
  HeatTransferSource<SourceFunction>::~HeatTransferSource()
  {
    return;
  }
  
  template< class SourceFunction >
  void HeatTransferSource<SourceFunction>::init_variables( libMesh::FEMSystem* system )
  {
     _T_var = system->add_variable( _T_var_name, this->_T_order, _T_FE_family);

    return;
  }

  template< class SourceFunction >
  void HeatTransferSource<SourceFunction>::element_time_derivative( bool /*compute_jacobian*/,
								    libMesh::FEMContext& context,
								    CachedValues& /*cache*/ )
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
    libMesh::DenseSubVector<libMesh::Number> &FT = *context.elem_subresiduals[_T_var]; // R_{T}

    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.element_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	libMesh::Real q = _source( x_qp[qp] );

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
