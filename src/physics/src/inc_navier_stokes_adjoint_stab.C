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
#include "grins/inc_navier_stokes_adjoint_stab.h"

//libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{

  IncompressibleNavierStokesAdjointStabilization::IncompressibleNavierStokesAdjointStabilization( const std::string& physics_name, 
												  const GetPot& input )
    : IncompressibleNavierStokesStabilizationBase(physics_name,input)
  {
    this->read_input_options(input);

    return;
  }

  IncompressibleNavierStokesAdjointStabilization::~IncompressibleNavierStokesAdjointStabilization()
  {
    return;
  }

  void IncompressibleNavierStokesAdjointStabilization::element_time_derivative( bool /*compute_jacobian*/,
										libMesh::FEMContext& context,
										CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("IncompressibleNavierStokesAdjointStabilization::element_time_derivative");
#endif

    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.dof_indices_var[this->_p_var].size();
    const unsigned int n_u_dofs = context.dof_indices_var[this->_u_var].size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.element_fe_var[this->_u_var]->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
      context.element_fe_var[this->_p_var]->get_dphi();

    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.element_fe_var[this->_u_var]->get_dphi();

    const std::vector<std::vector<libMesh::RealTensor> >& u_hessphi =
      context.element_fe_var[this->_u_var]->get_d2phi();

    libMesh::DenseSubVector<libMesh::Number> &Fu = *context.elem_subresiduals[this->_u_var]; // R_{p}
    libMesh::DenseSubVector<libMesh::Number> &Fv = *context.elem_subresiduals[this->_v_var]; // R_{p}
    libMesh::DenseSubVector<libMesh::Number> &Fw = *context.elem_subresiduals[this->_w_var]; // R_{w}
    libMesh::DenseSubVector<libMesh::Number> &Fp = *context.elem_subresiduals[this->_p_var]; // R_{p}

    unsigned int n_qpoints = context.element_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	libMesh::FEBase* fe = context.element_fe_var[this->_u_var];

	libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
	libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

	libMesh::RealGradient U( context.interior_value( this->_u_var, qp ),
				 context.interior_value( this->_v_var, qp ) );
	if( this->_dim == 3 )
	  U(2) = context.interior_value( this->_w_var, qp );
      
	libMesh::Real tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, this->_rho, U, this->_mu, this->_is_steady );
	libMesh::Real tau_C = this->_stab_helper.compute_tau_continuity( tau_M, g );

	libMesh::RealGradient RM_s = this->compute_res_momentum_steady( context, qp );
	libMesh::Real RC = compute_res_continuity( context, qp );

	// Now a loop over the pressure degrees of freedom.  This
	// computes the contributions of the continuity equation.
	for (unsigned int i=0; i != n_p_dofs; i++)
	  {
	    Fp(i) -= tau_M*RM_s*p_dphi[i][qp]*JxW[qp];
	  }

	for (unsigned int i=0; i != n_u_dofs; i++)
	  {
	    Fu(i) -= ( tau_M*RM_s(0)*this->_rho*U*u_gradphi[i][qp]
		       + tau_M*RM_s(0)*this->_mu*( u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) + u_hessphi[i][qp](2,2) )
		       + tau_C*RC*u_gradphi[i][qp](0) )*JxW[qp];

	    Fv(i) -= ( tau_M*RM_s(1)*this->_rho*U*u_gradphi[i][qp] 
		       + tau_M*RM_s(1)*this->_mu*( u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) + u_hessphi[i][qp](2,2) )
		       + tau_C*RC*u_gradphi[i][qp](1) )*JxW[qp];

	    if(this->_dim == 3)
	      Fw(i) -= ( tau_M*RM_s(2)*this->_rho*U*u_gradphi[i][qp] 
			 + tau_M*RM_s(2)*this->_mu*( u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) + u_hessphi[i][qp](2,2) )
			 + tau_C*RC*u_gradphi[i][qp](2) )*JxW[qp];
	  }

      }

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("IncompressibleNavierStokesAdjointStabilization::element_time_derivative");
#endif
    return;
  }

  void IncompressibleNavierStokesAdjointStabilization::mass_residual( bool /*compute_jacobian*/,
								      libMesh::FEMContext& context,
								      CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("IncompressibleNavierStokesAdjointStabilization::mass_residual");
#endif

    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.dof_indices_var[this->_p_var].size();
    const unsigned int n_u_dofs = context.dof_indices_var[this->_u_var].size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.element_fe_var[this->_u_var]->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
      context.element_fe_var[this->_p_var]->get_dphi();

    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.element_fe_var[this->_u_var]->get_dphi();

    const std::vector<std::vector<libMesh::RealTensor> >& u_hessphi =
      context.element_fe_var[this->_u_var]->get_d2phi();

    libMesh::DenseSubVector<libMesh::Number> &Fu = *context.elem_subresiduals[this->_u_var]; // R_{p}
    libMesh::DenseSubVector<libMesh::Number> &Fv = *context.elem_subresiduals[this->_v_var]; // R_{p}
    libMesh::DenseSubVector<libMesh::Number> &Fw = *context.elem_subresiduals[this->_w_var]; // R_{w}
    libMesh::DenseSubVector<libMesh::Number> &Fp = *context.elem_subresiduals[this->_p_var]; // R_{p}

    unsigned int n_qpoints = context.element_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	libMesh::FEBase* fe = context.element_fe_var[this->_u_var];

	libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
	libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

	libMesh::RealGradient U( context.fixed_interior_value( this->_u_var, qp ),
				 context.fixed_interior_value( this->_v_var, qp ) );
	if( this->_dim == 3 )
	  U(2) = context.fixed_interior_value( this->_w_var, qp );
      
	libMesh::Real tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, this->_rho, U, this->_mu, false );

	libMesh::RealGradient RM_t = this->compute_res_momentum_transient( context, qp );

	// Now a loop over the pressure degrees of freedom.  This
	// computes the contributions of the continuity equation.
	for (unsigned int i=0; i != n_p_dofs; i++)
	  {
	    Fp(i) += tau_M*RM_t*p_dphi[i][qp]*JxW[qp];
	  }

	for (unsigned int i=0; i != n_u_dofs; i++)
	  {
	    Fu(i) += tau_M*RM_t(0)*( this->_rho*U*u_gradphi[i][qp] 
				     + this->_mu*( u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) + u_hessphi[i][qp](2,2) ) 
				     )*JxW[qp];

	    Fv(i) += tau_M*RM_t(1)*( this->_rho*U*u_gradphi[i][qp]
				     + this->_mu*( u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) + u_hessphi[i][qp](2,2) ) 
				     )*JxW[qp];

	    Fw(i) += tau_M*RM_t(2)*( this->_rho*U*u_gradphi[i][qp]
				     + this->_mu*( u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) + u_hessphi[i][qp](2,2) ) 
				     )*JxW[qp];
	  }

      }

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("IncompressibleNavierStokesAdjointStabilization::mass_residual");
#endif
    return;
  }

} // namespace GRINS
