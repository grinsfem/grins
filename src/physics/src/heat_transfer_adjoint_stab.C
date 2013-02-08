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
#include "grins/heat_transfer_adjoint_stab.h"

// libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{

  HeatTransferAdjointStabilization::HeatTransferAdjointStabilization( const std::string& physics_name, 
								      const GetPot& input )
    : HeatTransferStabilizationBase(physics_name,input)
  {
    return;
  }

  HeatTransferAdjointStabilization::~HeatTransferAdjointStabilization()
  {
    return;
  }

  void HeatTransferAdjointStabilization::element_time_derivative( bool /*compute_jacobian*/,
								  libMesh::FEMContext& context,
								  CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("HeatTransferAdjointStabilization::element_time_derivative");
#endif

    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.dof_indices_var[this->_T_var].size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.element_fe_var[this->_T_var]->get_JxW();

    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.element_fe_var[this->_T_var]->get_dphi();

    const std::vector<std::vector<libMesh::RealTensor> >& T_hessphi =
      context.element_fe_var[this->_T_var]->get_d2phi();

    libMesh::DenseSubVector<libMesh::Number> &FT = *context.elem_subresiduals[this->_T_var]; // R_{T}

    unsigned int n_qpoints = context.element_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	libMesh::FEBase* fe = context.element_fe_var[this->_T_var];

	libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
	libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

	libMesh::RealGradient U( context.interior_value( this->_u_var, qp ),
				 context.interior_value( this->_v_var, qp ) );
	if( this->_dim == 3 )
	  U(2) = context.interior_value( this->_w_var, qp );
      
	libMesh::Real tau_E = this->_stab_helper.compute_tau_energy( context, G, _rho, _Cp, _k,  U, this->_is_steady );

	libMesh::Real RE_s = this->compute_res_steady( context, qp );

	for (unsigned int i=0; i != n_T_dofs; i++)
	  {
	    FT(i) += tau_E*RE_s*( _rho*_Cp*U*T_gradphi[i][qp]
				  + _k*(T_hessphi[i][qp](0,0) + T_hessphi[i][qp](1,1) + T_hessphi[i][qp](2,2)) 
				  )*JxW[qp];
	  }

      }

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("HeatTransferAdjointStabilization::element_time_derivative");
#endif
    return;
  }

  void HeatTransferAdjointStabilization::mass_residual( bool /*compute_jacobian*/,
							libMesh::FEMContext& context,
							CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("HeatTransferAdjointStabilization::mass_residual");
#endif

    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.dof_indices_var[this->_T_var].size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.element_fe_var[this->_T_var]->get_JxW();

    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.element_fe_var[this->_T_var]->get_dphi();

    const std::vector<std::vector<libMesh::RealTensor> >& T_hessphi =
      context.element_fe_var[this->_T_var]->get_d2phi();

    libMesh::DenseSubVector<libMesh::Number> &FT = *context.elem_subresiduals[this->_T_var]; // R_{T}

    unsigned int n_qpoints = context.element_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	libMesh::FEBase* fe = context.element_fe_var[this->_T_var];

	libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
	libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

	libMesh::RealGradient U( context.fixed_interior_value( this->_u_var, qp ),
				 context.fixed_interior_value( this->_v_var, qp ) );
	if( this->_dim == 3 )
	  U(2) = context.fixed_interior_value( this->_w_var, qp );
      
	libMesh::Real tau_E = this->_stab_helper.compute_tau_energy( context, G, _rho, _Cp, _k,  U, false );

	libMesh::Real RE_t = this->compute_res_transient( context, qp );

	for (unsigned int i=0; i != n_T_dofs; i++)
	  {
	    FT(i) -= tau_E*RE_t*( _rho*_Cp*U*T_gradphi[i][qp]
				  + _k*(T_hessphi[i][qp](0,0) + T_hessphi[i][qp](1,1) + T_hessphi[i][qp](2,2)) 
				  )*JxW[qp];
	  }

      }

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("HeatTransferAdjointStabilization::mass_residual");
#endif
    return;
  }

} // namespace GRINS
