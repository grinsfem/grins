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

#include "inc_navier_stokes_adjoint_stab.h"

GRINS::IncompressibleNavierStokesAdjointStabilization::IncompressibleNavierStokesAdjointStabilization( const std::string& physics_name, 
											   const GetPot& input )
  : GRINS::IncompressibleNavierStokesStabilizationBase(physics_name,input)
{
  this->read_input_options(input);

  return;
}

GRINS::IncompressibleNavierStokesAdjointStabilization::~IncompressibleNavierStokesAdjointStabilization()
{
  return;
}

void GRINS::IncompressibleNavierStokesAdjointStabilization::read_input_options( const GetPot& input )
{
  return;
}

bool GRINS::IncompressibleNavierStokesAdjointStabilization::element_time_derivative( bool request_jacobian,
										     libMesh::DiffContext& context,
										     libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("IncompressibleNavierStokesAdjointStabilization::element_time_derivative");
#endif

  libMesh::FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // The number of local degrees of freedom in each variable.
  const unsigned int n_p_dofs = c.dof_indices_var[this->_p_var].size();
  const unsigned int n_u_dofs = c.dof_indices_var[this->_u_var].size();

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[this->_u_var]->get_JxW();

  // The pressure shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
    c.element_fe_var[this->_p_var]->get_dphi();

  const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
    c.element_fe_var[this->_u_var]->get_dphi();

  const std::vector<std::vector<libMesh::RealTensor> >& u_hessphi =
    c.element_fe_var[this->_u_var]->get_d2phi();

  libMesh::DenseSubVector<Number> &Fu = *c.elem_subresiduals[this->_u_var]; // R_{p}
  libMesh::DenseSubVector<Number> &Fv = *c.elem_subresiduals[this->_v_var]; // R_{p}
  libMesh::DenseSubVector<Number> &Fp = *c.elem_subresiduals[this->_p_var]; // R_{p}

  unsigned int n_qpoints = c.element_qrule->n_points();

  bool is_steady = (system->time_solver)->is_steady();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      libMesh::FEBase* fe = c.element_fe_var[this->_u_var];

      libMesh::RealGradient g = this->_stab_helper.compute_g( fe, c, qp );
      libMesh::RealTensor G = this->_stab_helper.compute_G( fe, c, qp );

      libMesh::RealGradient U( c.interior_value( this->_u_var, qp ),
			       c.interior_value( this->_v_var, qp ) );
      if( this->_dim == 3 )
	U(2) = c.interior_value( this->_w_var, qp );
      
      libMesh::Real tau_M = this->_stab_helper.compute_tau_momentum( c, qp, g, G, this->_rho, U, this->_mu, is_steady );
      libMesh::Real tau_C = this->_stab_helper.compute_tau_continuity( tau_M, g );

      libMesh::RealGradient RM_s = this->compute_res_momentum_steady( c, qp );
      libMesh::Real RC = compute_res_continuity( c, qp );

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) -= tau_M*RM_s*p_dphi[i][qp]*JxW[qp];
	}

      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fu(i) -= ( tau_M*RM_s(0)*this->_rho*U*u_gradphi[i][qp]
		     + this->_mu*( u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) )
		     + tau_C*RC*u_gradphi[i][qp](0) )*JxW[qp];

	  Fv(i) -= ( tau_M*RM_s(1)*this->_rho*U*u_gradphi[i][qp] 
		     + this->_mu*( u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) )
		     + tau_C*RC*u_gradphi[i][qp](1) )*JxW[qp];
	}

    }

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("IncompressibleNavierStokesAdjointStabilization::element_time_derivative");
#endif
  return request_jacobian;
}

bool GRINS::IncompressibleNavierStokesAdjointStabilization::mass_residual( bool request_jacobian,
									   libMesh::DiffContext& context,
									   libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("IncompressibleNavierStokesAdjointStabilization::mass_residual");
#endif

  libMesh::FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // The number of local degrees of freedom in each variable.
  const unsigned int n_p_dofs = c.dof_indices_var[this->_p_var].size();
  const unsigned int n_u_dofs = c.dof_indices_var[this->_u_var].size();

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[this->_u_var]->get_JxW();

  // The pressure shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
    c.element_fe_var[this->_p_var]->get_dphi();

  const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
    c.element_fe_var[this->_u_var]->get_dphi();

  const std::vector<std::vector<libMesh::RealTensor> >& u_hessphi =
    c.element_fe_var[this->_u_var]->get_d2phi();

  libMesh::DenseSubVector<Number> &Fu = *c.elem_subresiduals[this->_u_var]; // R_{p}
  libMesh::DenseSubVector<Number> &Fv = *c.elem_subresiduals[this->_v_var]; // R_{p}
  libMesh::DenseSubVector<Number> &Fp = *c.elem_subresiduals[this->_p_var]; // R_{p}

  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      libMesh::FEBase* fe = c.element_fe_var[this->_u_var];

      libMesh::RealGradient g = this->_stab_helper.compute_g( fe, c, qp );
      libMesh::RealTensor G = this->_stab_helper.compute_G( fe, c, qp );

      libMesh::RealGradient U( c.fixed_interior_value( this->_u_var, qp ),
			       c.fixed_interior_value( this->_v_var, qp ) );
      if( this->_dim == 3 )
	U(2) = c.fixed_interior_value( this->_w_var, qp );
      
      libMesh::Real tau_M = this->_stab_helper.compute_tau_momentum( c, qp, g, G, this->_rho, U, this->_mu, false );
      libMesh::Real tau_C = this->_stab_helper.compute_tau_continuity( tau_M, g );

      libMesh::RealGradient RM_t = this->compute_res_momentum_transient( c, qp );

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) += tau_M*RM_t*p_dphi[i][qp]*JxW[qp];
	}

      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fu(i) += tau_M*RM_t(0)*( this->_rho*U*u_gradphi[i][qp] 
				   + this->_mu*( u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) ) 
				   )*JxW[qp];

	  Fv(i) += tau_M*RM_t(1)*( this->_rho*U*u_gradphi[i][qp]
				   + this->_mu*( u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) ) 
				   )*JxW[qp];
	}

    }

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("IncompressibleNavierStokesAdjointStabilization::mass_residual");
#endif
  return request_jacobian;
}

bool GRINS::IncompressibleNavierStokesAdjointStabilization::element_constraint( bool request_jacobian,
										libMesh::DiffContext& context,
										libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  //this->_timer->BeginTimer("IncompressibleNavierStokesAdjointStabilization::element_constraint");
#endif

  //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

#ifdef USE_GRVY_TIMERS
  //this->_timer->EndTimer("IncompressibleNavierStokesAdjointStabilization::element_constraint");
#endif

  return request_jacobian;
}

bool GRINS::IncompressibleNavierStokesAdjointStabilization::side_time_derivative( bool request_jacobian,
										  libMesh::DiffContext& context,
										  libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
      //this->_timer->BeginTimer("IncompressibleNavierStokesAdjointStabilization::side_time_derivative");
#endif
      //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

#ifdef USE_GRVY_TIMERS
      //this->_timer->EndTimer("IncompressibleNavierStokesAdjointStabilization::side_time_derivative");
#endif

  return request_jacobian;
}

bool GRINS::IncompressibleNavierStokesAdjointStabilization::side_constraint( bool request_jacobian,
									     libMesh::DiffContext& context,
									     libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  //this->_timer->BeginTimer("IncompressibleNavierStokesAdjointStabilization::side_constraint");
#endif

  //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

#ifdef USE_GRVY_TIMERS
  //this->_timer->EndTimer("IncompressibleNavierStokesAdjointStabilization::side_constraint");
#endif

  return request_jacobian;
}
