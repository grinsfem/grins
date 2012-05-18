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

#include "low_mach_navier_stokes_vms_stab.h"

template<class Mu, class SH, class TC>
GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::LowMachNavierStokesVMSStabilization( const std::string& physics_name, 
											   const GetPot& input )
  : GRINS::LowMachNavierStokesBase<Mu,SH,TC>(physics_name,input)
{
  this->read_input_options(input);

  return;
}

template<class Mu, class SH, class TC>
GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::~LowMachNavierStokesVMSStabilization()
{
  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::read_input_options( const GetPot& input )
{
  this->_C = input("Physics/"+this->_physics_name+"/tau_constant", 1 );
  this->_tau_factor = input("Physics/"+this->_physics_name+"/tau_factor", 0.5 );
  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::init_context( libMesh::DiffContext &context )
{
  // First call base class
  GRINS::LowMachNavierStokesBase<Mu,SH,TC>::init_context(context);

  libMesh::FEMContext &c = libmesh_cast_ref<libMesh::FEMContext&>(context);
  
  // We need pressure derivatives
  c.element_fe_var[this->_p_var]->get_dphi();

  // We also need second derivatives, so initialize those.
  c.element_fe_var[this->_u_var]->get_d2phi();
  c.element_fe_var[this->_T_var]->get_d2phi();

  return;
}

template<class Mu, class SH, class TC>
bool GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::element_time_derivative( bool request_jacobian,
										    libMesh::DiffContext& context,
										    libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("LowMachNavierStokesVMSStabilization::element_time_derivative");
#endif

  libMesh::FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  this->assemble_continuity_time_deriv( request_jacobian, c, system );
  this->assemble_momentum_time_deriv( request_jacobian, c, system );
  this->assemble_energy_time_deriv( request_jacobian, c, system );

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("LowMachNavierStokesVMSStabilization::element_time_derivative");
#endif
  return request_jacobian;
}

template<class Mu, class SH, class TC>
bool GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::mass_residual( bool request_jacobian,
									  libMesh::DiffContext& context,
									  libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("LowMachNavierStokesVMSStabilization::mass_residual");
#endif

  libMesh::FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  this->assemble_continuity_mass_residual( request_jacobian, c, system );
  this->assemble_momentum_mass_residual( request_jacobian, c, system );
  this->assemble_energy_mass_residual( request_jacobian, c, system );

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("LowMachNavierStokesVMSStabilization::mass_residual");
#endif
  return request_jacobian;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::assemble_continuity_time_deriv( bool request_jacobian,
											   libMesh::FEMContext& c,
											   libMesh::FEMSystem* system )
{
  // The number of local degrees of freedom in each variable.
  const unsigned int n_p_dofs = c.dof_indices_var[this->_p_var].size();

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[this->_u_var]->get_JxW();

  // The pressure shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
    c.element_fe_var[this->_p_var]->get_dphi();

  libMesh::DenseSubVector<Number> &Fp = *c.elem_subresiduals[this->_p_var]; // R_{p}

  unsigned int n_qpoints = c.element_qrule->n_points();

  bool is_steady = (system->time_solver)->is_steady();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      libMesh::RealGradient g = this->compute_g( c, qp );
      libMesh::RealTensor G = this->compute_G( c, qp );

      libMesh::Real T = c.interior_value( this->_T_var, qp );
      libMesh::Real rho = this->compute_rho( T, this->get_p0_steady( c, qp ) );

      libMesh::RealGradient U( c.interior_value( this->_u_var, qp ),
			       c.interior_value( this->_v_var, qp ) );
      if( this->_dim == 3 )
	U(2) = c.interior_value( this->_w_var, qp );

      libMesh::Real tau_M = this->compute_tau_momentum( c, qp, g, G, rho, U, T, is_steady );
      libMesh::Real tau_E = this->compute_tau_energy( c, qp, g, G, rho, U, T, is_steady );

      libMesh::RealGradient RM_s = this->compute_res_momentum_steady( c, qp );
      libMesh::Real RE_s = this->compute_res_energy_steady( c, qp );

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) += tau_M*RM_s*p_dphi[i][qp]*JxW[qp];
	  Fp(i) -= tau_E*RE_s*(U*p_dphi[i][qp])/T*JxW[qp];
	}

    }

  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::assemble_momentum_time_deriv( bool request_jacobian,
											 libMesh::FEMContext& c,
											 libMesh::FEMSystem* system )
{
  // The number of local degrees of freedom in each variable.
  const unsigned int n_u_dofs = c.dof_indices_var[this->_u_var].size();

  // Check number of dofs is same for _u_var, v_var and w_var.
  libmesh_assert (n_u_dofs == c.dof_indices_var[this->_v_var].size());
  if (this->_dim == 3)
    libmesh_assert (n_u_dofs == c.dof_indices_var[this->_w_var].size());

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[this->_u_var]->get_JxW();

  // The pressure shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& u_phi =
    c.element_fe_var[this->_u_var]->get_phi();

// The velocity shape function gradients at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
    c.element_fe_var[this->_u_var]->get_dphi();

  const std::vector<std::vector<libMesh::RealTensor> >& u_hessphi =
    c.element_fe_var[this->_u_var]->get_d2phi();

  libMesh::DenseSubVector<Number> &Fu = *c.elem_subresiduals[this->_u_var]; // R_{u}
  libMesh::DenseSubVector<Number> &Fv = *c.elem_subresiduals[this->_v_var]; // R_{v}
  libMesh::DenseSubVector<Number> &Fw = *c.elem_subresiduals[this->_w_var]; // R_{w}

  unsigned int n_qpoints = c.element_qrule->n_points();

  bool is_steady = (system->time_solver)->is_steady();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      libMesh::Real T = c.interior_value( this->_T_var, qp );
      libMesh::Real rho = this->compute_rho( T, this->get_p0_steady( c, qp ) );

      libMesh::RealGradient U( c.interior_value(this->_u_var, qp),
			       c.interior_value(this->_v_var, qp) );

      libMesh::RealGradient grad_u = c.interior_gradient(this->_u_var, qp);
      libMesh::RealGradient grad_v = c.interior_gradient(this->_v_var, qp);
      libMesh::RealGradient grad_w;

      if( this->_dim == 3 )
	{
	  U(2) = c.interior_value(this->_w_var, qp);
	  grad_w = c.interior_gradient(this->_w_var, qp);
	}

      libMesh::RealGradient g = this->compute_g( c, qp );
      libMesh::RealTensor G = this->compute_G( c, qp );

      libMesh::Real tau_M = this->compute_tau_momentum( c, qp, g, G, rho, U, T, is_steady );
      libMesh::Real tau_C = this->compute_tau_continuity( tau_M, g, G, U, rho );

      libMesh::Real RC_s = this->compute_res_continuity_steady( c, qp );
      libMesh::RealGradient RM_s = this->compute_res_momentum_steady( c, qp );

      /*
      std::cout << "g = " << g << std::endl
		<< "G = " << G << std::endl;
      */

      /*
      std::cout << "tau_M = " << tau_M << ", tau_C = " << tau_C << std::endl
		<< "RC_s = " << RC_s << std::endl
		<< "RM_s = " << RM_s << std::endl;
      */
      libMesh::Real mu = this->_mu(T);

      for (unsigned int i=0; i != n_u_dofs; i++)
        {
	  Fu(i) += ( tau_C*RC_s*u_gradphi[i][qp](0)
		     //+ rho*tau_M*RM_s*grad_u*u_phi[i][qp]
		     + tau_M*RM_s(0)*rho*U*u_gradphi[i][qp] 
		     + mu*tau_M*RM_s(0)*(u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) 
					 + u_hessphi[i][qp](0,0) + u_hessphi[i][qp](0,1) - 2.0/3.0*(u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,0)) 
					 ) )*JxW[qp];
		     //+ tau_M*RM_s(0)*rho*tau_M*RM_s*u_gradphi[i][qp] )*JxW[qp];

	  Fv(i) += ( tau_C*RC_s*u_gradphi[i][qp](1)
		     //+ rho*tau_M*RM_s*grad_v*u_phi[i][qp]
		     + tau_M*RM_s(1)*rho*U*u_gradphi[i][qp]
		     + mu*tau_M*RM_s(1)*(u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) 
					 + u_hessphi[i][qp](1,0) + u_hessphi[i][qp](1,1) - 2.0/3.0*(u_hessphi[i][qp](0,1) + u_hessphi[i][qp](1,1)) 
					 ) )*JxW[qp];
	  //+ tau_M*RM_s(1)*rho*tau_M*RM_s*u_gradphi[i][qp] )*JxW[qp];

	  if( this->_dim == 3 )
	    {
	      Fw(i) += ( -tau_C*RC_s*u_gradphi[i][qp](2)
			 + rho*tau_M*RM_s*grad_w*u_phi[i][qp]
			 - tau_M*RM_s(2)*rho*U*u_gradphi[i][qp]
			 + tau_M*RM_s(2)*rho*tau_M*RM_s*u_gradphi[i][qp] )*JxW[qp];
	    }
	}

    }
  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::assemble_energy_time_deriv( bool request_jacobian,
										       libMesh::FEMContext& c,
										       libMesh::FEMSystem* system )
{
  // The number of local degrees of freedom in each variable.
  const unsigned int n_T_dofs = c.dof_indices_var[this->_T_var].size();

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[this->_T_var]->get_JxW();

  // The temperature shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& T_phi =
    c.element_fe_var[this->_T_var]->get_phi();

  // The temperature shape functions gradients at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
    c.element_fe_var[this->_T_var]->get_dphi();

  const std::vector<std::vector<libMesh::RealTensor> >& T_hessphi =
    c.element_fe_var[this->_T_var]->get_d2phi();

  libMesh::DenseSubVector<Number> &FT = *c.elem_subresiduals[this->_T_var]; // R_{T}

  unsigned int n_qpoints = c.element_qrule->n_points();

  bool is_steady = (system->time_solver)->is_steady();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      libMesh::Number u, v, w;
      u = c.interior_value(this->_u_var, qp);
      v = c.interior_value(this->_v_var, qp);
      if (this->_dim == 3)
        w = c.interior_value(this->_w_var, qp);

      libMesh::Gradient grad_T = c.interior_gradient(this->_T_var, qp);

      libMesh::NumberVectorValue U(u,v);
      if (this->_dim == 3)
        U(2) = w;
      
      libMesh::Real T = c.interior_value( this->_T_var, qp );
      libMesh::Real rho = this->compute_rho( T, this->get_p0_steady( c, qp ) );

      libMesh::Number rho_cp = rho*this->_cp(T);

      libMesh::RealGradient g = this->compute_g( c, qp );
      libMesh::RealTensor G = this->compute_G( c, qp );

      libMesh::Real tau_M = this->compute_tau_momentum( c, qp, g, G, rho, U, T, is_steady );
      libMesh::Real tau_E = this->compute_tau_energy( c, qp, g, G, rho, U, T, is_steady );

      libMesh::Real RE_s = this->compute_res_energy_steady( c, qp );
      libMesh::RealGradient RM_s = this->compute_res_momentum_steady( c, qp );

      //std::cout << "tau_E = " << tau_E << ", RE_s = " << RE_s << std::endl;

      libMesh::Real k = this->_k(T);

      for (unsigned int i=0; i != n_T_dofs; i++)
        {
          FT(i) += ( //rho_cp*tau_M*RM_s*grad_T*T_phi[i][qp] 
		    + rho_cp*tau_E*RE_s*U*T_gradphi[i][qp]
		    + tau_E*RE_s*k*(T_hessphi[i][qp](0,0) + T_hessphi[i][qp](1,1) + T_hessphi[i][qp](2,2) ) )*JxW[qp];
	    //+ rho_cp*tau_E*RE_s*tau_M*RM_s*T_gradphi[i][qp] )*JxW[qp];
	}

    }

  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::assemble_continuity_mass_residual( bool request_jacobian,
											      libMesh::FEMContext& c,
											      libMesh::FEMSystem* system )
{
  // The number of local degrees of freedom in each variable.
  const unsigned int n_p_dofs = c.dof_indices_var[this->_p_var].size();

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[this->_u_var]->get_JxW();

  // The pressure shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
    c.element_fe_var[this->_p_var]->get_dphi();

  libMesh::DenseSubVector<Number> &Fp = *c.elem_subresiduals[this->_p_var]; // R_{p}

  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      libMesh::RealGradient g = this->compute_g( c, qp );
      libMesh::RealTensor G = this->compute_G( c, qp );

      libMesh::Real T = c.fixed_interior_value( this->_T_var, qp );
      libMesh::Real rho = this->compute_rho( T, this->get_p0_transient( c, qp ) );

      libMesh::RealGradient U( c.fixed_interior_value( this->_u_var, qp ),
			       c.fixed_interior_value( this->_v_var, qp ) );
      if( this->_dim == 3 )
	U(2) = c.fixed_interior_value( this->_w_var, qp );

      libMesh::Real tau_M = this->compute_tau_momentum( c, qp, g, G, rho, U, T, false );
      libMesh::RealGradient RM_t = this->compute_res_momentum_transient( c, qp );

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) -= tau_M*RM_t*p_dphi[i][qp]*JxW[qp];
	}
    }

  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::assemble_momentum_mass_residual( bool request_jacobian,
											    libMesh::FEMContext& c,
											    libMesh::FEMSystem* system )
{
  // The number of local degrees of freedom in each variable.
  const unsigned int n_u_dofs = c.dof_indices_var[this->_u_var].size();

  // Check number of dofs is same for _u_var, v_var and w_var.
  libmesh_assert (n_u_dofs == c.dof_indices_var[this->_v_var].size());
  if (this->_dim == 3)
    libmesh_assert (n_u_dofs == c.dof_indices_var[this->_w_var].size());

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[this->_u_var]->get_JxW();

  // The pressure shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& u_phi =
    c.element_fe_var[this->_u_var]->get_phi();

// The velocity shape function gradients at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
    c.element_fe_var[this->_u_var]->get_dphi();

  libMesh::DenseSubVector<Number> &Fu = *c.elem_subresiduals[this->_u_var]; // R_{u}
  libMesh::DenseSubVector<Number> &Fv = *c.elem_subresiduals[this->_v_var]; // R_{v}
  libMesh::DenseSubVector<Number> &Fw = *c.elem_subresiduals[this->_w_var]; // R_{w}

  unsigned int n_qpoints = c.element_qrule->n_points();
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      libMesh::Real T = c.fixed_interior_value( this->_T_var, qp );
      libMesh::Real rho = this->compute_rho( T, this->get_p0_transient( c, qp ) );

      libMesh::RealGradient U( c.fixed_interior_value(this->_u_var, qp),
			       c.fixed_interior_value(this->_v_var, qp) );

      libMesh::RealGradient grad_u = c.fixed_interior_gradient(this->_u_var, qp);
      libMesh::RealGradient grad_v = c.fixed_interior_gradient(this->_v_var, qp);
      libMesh::RealGradient grad_w;

      if( this->_dim == 3 )
	{
	  U(2) = c.fixed_interior_value(this->_w_var, qp);
	  grad_w = c.fixed_interior_gradient(this->_w_var, qp);
	}

      libMesh::RealGradient g = this->compute_g( c, qp );
      libMesh::RealTensor G = this->compute_G( c, qp );

      libMesh::Real tau_M = this->compute_tau_momentum( c, qp, g, G, rho, U, T, false );
      libMesh::Real tau_C = this->compute_tau_continuity( tau_M, g, G, U, rho );

      libMesh::Real RC_t = this->compute_res_continuity_transient( c, qp );
      libMesh::RealGradient RM_s = this->compute_res_momentum_steady( c, qp );
      libMesh::RealGradient RM_t = this->compute_res_momentum_transient( c, qp );
      
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
	  Fu(i) += ( tau_C*RC_t*u_gradphi[i][qp](0)
		     - rho*tau_M*RM_t*grad_u*u_phi[i][qp]
		     + tau_M*RM_t(0)*rho*U*u_gradphi[i][qp]
		     - tau_M*(RM_s(0)+RM_t(0))*rho*tau_M*RM_t*u_gradphi[i][qp]
		     - tau_M*RM_t(0)*rho*tau_M*RM_s*u_gradphi[i][qp] )*JxW[qp];

	  Fv(i) += ( tau_C*RC_t*u_gradphi[i][qp](1)
		     - rho*tau_M*RM_t*grad_v*u_phi[i][qp]
		     + tau_M*RM_t(1)*rho*U*u_gradphi[i][qp]
		     - tau_M*(RM_s(1)+RM_t(1))*rho*tau_M*RM_t*u_gradphi[i][qp]
		     - tau_M*RM_t(1)*rho*tau_M*RM_s*u_gradphi[i][qp] )*JxW[qp];

	  if( this->_dim == 3 )
	    {
	      Fw(i) += ( tau_C*RC_t*u_gradphi[i][qp](2)
			 - rho*tau_M*RM_t*grad_w*u_phi[i][qp]
			 + tau_M*RM_t(2)*rho*U*u_gradphi[i][qp]
			 - tau_M*(RM_s(2)+RM_t(2))*rho*tau_M*RM_t*u_gradphi[i][qp]
			 - tau_M*RM_t(2)*rho*tau_M*RM_s*u_gradphi[i][qp] )*JxW[qp];
	    }
	}

    }
  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::assemble_energy_mass_residual( bool request_jacobian,
											  libMesh::FEMContext& c,
											  libMesh::FEMSystem* system )
{
  // The number of local degrees of freedom in each variable.
  const unsigned int n_T_dofs = c.dof_indices_var[this->_T_var].size();

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[this->_T_var]->get_JxW();

  // The temperature shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& T_phi =
    c.element_fe_var[this->_T_var]->get_phi();

  // The temperature shape functions gradients at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
    c.element_fe_var[this->_T_var]->get_dphi();

  libMesh::DenseSubVector<Number> &FT = *c.elem_subresiduals[this->_T_var]; // R_{T}

  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      libMesh::Number u, v, w;
      u = c.fixed_interior_value(this->_u_var, qp);
      v = c.fixed_interior_value(this->_v_var, qp);
      if (this->_dim == 3)
        w = c.fixed_interior_value(this->_w_var, qp);

      libMesh::Gradient grad_T = c.fixed_interior_gradient(this->_T_var, qp);

      libMesh::NumberVectorValue U(u,v);
      if (this->_dim == 3)
        U(2) = w;

      libMesh::Real T = c.fixed_interior_value( this->_T_var, qp );
      libMesh::Real rho = this->compute_rho( T, this->get_p0_transient( c, qp ) );

      libMesh::Number rho_cp = rho*this->_cp(T);

      libMesh::RealGradient g = this->compute_g( c, qp );
      libMesh::RealTensor G = this->compute_G( c, qp );

      libMesh::Real tau_M = this->compute_tau_momentum( c, qp, g, G, rho, U, T, false );
      libMesh::Real tau_E = this->compute_tau_energy( c, qp, g, G, rho, U, T, false );

      libMesh::Real RE_s = this->compute_res_energy_steady( c, qp );
      libMesh::Real RE_t = this->compute_res_energy_transient( c, qp );

      libMesh::RealGradient RM_s = this->compute_res_momentum_steady( c, qp );
      libMesh::RealGradient RM_t = this->compute_res_momentum_transient( c, qp );

      for (unsigned int i=0; i != n_T_dofs; i++)
        {
          FT(i) += ( -rho_cp*tau_M*RM_t*grad_T*T_phi[i][qp] 
		     +rho_cp*tau_E*RE_t*U*T_gradphi[i][qp]
		     - rho_cp*tau_E*(RE_s+RE_t)*tau_M*RM_t*T_gradphi[i][qp]
		     - rho_cp*tau_E*RE_t*tau_M*RM_s*T_gradphi[i][qp] )*JxW[qp];
	}

    }

  return;
}

template<class Mu, class SH, class TC>
libMesh::Real GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::compute_res_continuity_steady( libMesh::FEMContext& c,
												   unsigned int qp ) const
{
  libMesh::Real T = c.interior_value(this->_T_var, qp);
  libMesh::RealGradient grad_T = c.interior_gradient(this->_T_var, qp);

  libMesh::RealGradient U( c.interior_value(this->_u_var, qp),
			   c.interior_value(this->_v_var, qp) );  

  libMesh::RealGradient grad_u, grad_v;

  grad_u = c.interior_gradient(this->_u_var, qp);
  grad_v = c.interior_gradient(this->_v_var, qp);

  libMesh::Real divU = grad_u(0) + grad_v(1);

  if( this->_dim == 3 )
    {
      U(2) = c.interior_value(this->_w_var, qp);
      divU += (c.interior_gradient(this->_w_var, qp))(2);
    }

  return divU - (U*grad_T)/T;
}

template<class Mu, class SH, class TC>
libMesh::Real GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::compute_res_continuity_transient( libMesh::FEMContext& c,
												      unsigned int qp ) const
{
  libMesh::Real T = c.fixed_interior_value(this->_T_var, qp);
  libMesh::Real T_dot = c.interior_value(this->_T_var, qp);

  libMesh::Real RC_t = -T_dot/T;

  if( this->_enable_thermo_press_calc )
    {
      libMesh::Real p0 = c.fixed_interior_value(this->_p0_var, qp);
      libMesh::Real p0_dot = c.interior_value(this->_p0_var, qp);

      RC_t += p0_dot/p0;
    }

  return RC_t;
}

template<class Mu, class SH, class TC>
libMesh::RealGradient GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::compute_res_momentum_steady( libMesh::FEMContext& c,
													 unsigned int qp ) const
{
  libMesh::Real T = c.interior_value(this->_T_var, qp);

  libMesh::Real rho = this->compute_rho(T, this->get_p0_steady(c,qp) );

  libMesh::RealGradient rhoU( rho*c.interior_value(this->_u_var, qp), 
			      rho*c.interior_value(this->_v_var, qp) );
  if(this->_dim == 3)
    rhoU(2) = rho*c.interior_value(this->_w_var, qp);

  libMesh::RealGradient rhoUdotGradU( rhoU*c.interior_gradient(this->_u_var, qp),
				      rhoU*c.interior_gradient(this->_v_var, qp) );
  if(this->_dim == 3)
    rhoUdotGradU(2) = rhoU*c.interior_gradient(this->_w_var, qp);

  libMesh::RealGradient grad_p = c.interior_gradient(this->_p_var, qp);

  libMesh::RealTensor hess_u = c.interior_hessian(this->_u_var, qp);
  libMesh::RealTensor hess_v = c.interior_hessian(this->_v_var, qp);

  libMesh::RealGradient divGradU( hess_u(0,0) + hess_u(1,1),
				  hess_v(0,0) + hess_v(1,1) );

  libMesh::RealGradient divGradUT( hess_u(0,0) + hess_v(0,1),
				   hess_u(1,0) + hess_v(1,1) );

  libMesh::RealGradient divdivU( hess_u(0,0) + hess_v(1,0),
				 hess_u(0,1) + hess_v(1,1) );

  if(this->_dim == 3)
    {
      libMesh::RealTensor hess_w = c.interior_hessian(this->_w_var, qp);

      divGradU(0) += hess_u(2,2);
      divGradU(1) += hess_v(2,2);
      divGradU(2) = hess_w(0,0) + hess_w(1,1) + hess_w(2,2);

      divGradUT(0) += hess_w(0,2);
      divGradUT(1) += hess_w(1,2);
      divGradUT(2) = hess_u(2,0) + hess_v(2,1) + hess_w(2,2);

      divdivU(0) += hess_w(2,0);
      divdivU(1) += hess_w(2,1);
      divdivU(2) = hess_u(0,2) + hess_v(1,2) + hess_w(2,2);
    }

  libMesh::RealGradient divT = this->_mu(T)*(divGradU + divGradUT - 2.0/3.0*divdivU);

  if( this->_mu.deriv(T) != 0.0 )
    {
      libMesh::Gradient grad_T = c.interior_gradient(this->_T_var, qp);

      libMesh::Gradient grad_u = c.interior_gradient(this->_u_var, qp);
      libMesh::Gradient grad_v = c.interior_gradient(this->_v_var, qp);

      libMesh::Gradient gradTgradu( grad_T*grad_u, grad_T*grad_v );

      libMesh::Gradient gradTgraduT( grad_T(0)*grad_u(0) + grad_T(1)*grad_u(1),
				     grad_T(0)*grad_v(0) + grad_T(1)*grad_v(1) );

      libMesh::Real divU = grad_u(0) + grad_v(1);

      libMesh::Gradient gradTdivU( grad_T(0)*divU, grad_T(1)*divU );

      if(this->_dim == 3)
	{
	  libMesh::Gradient grad_w = c.interior_gradient(this->_w_var, qp);

	  gradTgradu(2) = grad_T*grad_w;

	  gradTgraduT(0) += grad_T(2)*grad_u(2);
	  gradTgraduT(1) += grad_T(2)*grad_v(2);
	  gradTgraduT(2) = grad_T(0)*grad_w(0) + grad_T(1)*grad_w(1) + grad_T(2)*grad_w(2);

	  divU += grad_w(2);
	  gradTdivU(0) += grad_T(0)*grad_w(2);
	  gradTdivU(1) += grad_T(1)*grad_w(2);
	  gradTdivU(2) += grad_T(2)*divU;
	}
      
      divT += this->_mu.deriv(T)*( gradTgradu + gradTgraduT - 2.0/3.0*gradTdivU );
    }

  libMesh::RealGradient rhog( rho*this->_g(0), rho*this->_g(1) );
  if(this->_dim == 3)
    rhog(2) = rho*this->_g(2);

  return rhoUdotGradU + grad_p - divT - rhog;
}

template<class Mu, class SH, class TC>
libMesh::RealGradient GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::compute_res_momentum_transient( libMesh::FEMContext& c,
													    unsigned int qp ) const
{
  libMesh::Real T = c.fixed_interior_value(this->_T_var, qp);
  libMesh::Real rho = this->compute_rho(T, this->get_p0_transient(c,qp) );

  libMesh::RealGradient u_dot( c.interior_value(this->_u_var, qp), c.interior_value(this->_v_var, qp) );

  if(this->_dim == 3)
    u_dot(2) = c.interior_value(this->_w_var, qp);

  return rho*u_dot;
}

template<class Mu, class SH, class TC>
libMesh::Real GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::compute_res_energy_steady( libMesh::FEMContext& c,
											       unsigned int qp ) const
{
  libMesh::Real T = c.interior_value(this->_T_var, qp);
  libMesh::Gradient grad_T = c.interior_gradient(this->_T_var, qp);
  libMesh::Tensor hess_T = c.interior_hessian(this->_T_var, qp);

  libMesh::Real rho = this->compute_rho(T, this->get_p0_steady(c,qp) );
  libMesh::Real rho_cp = rho*this->_cp(T);

  libMesh::RealGradient rhocpU( rho_cp*c.interior_value(this->_u_var, qp), 
				rho_cp*c.interior_value(this->_v_var, qp) );
  if(this->_dim == 3)
    rhocpU(2) = rho_cp*c.interior_value(this->_w_var, qp);

  return rhocpU*grad_T - this->_k.deriv(T)*(grad_T*grad_T) - this->_k(T)*(hess_T(0,0) + hess_T(1,1) + hess_T(2,2));
}


template<class Mu, class SH, class TC>
libMesh::Real GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::compute_res_energy_transient( libMesh::FEMContext& c,
												  unsigned int qp ) const
{
  libMesh::Real T = c.fixed_interior_value(this->_T_var, qp);
  libMesh::Real rho = this->compute_rho(T, this->get_p0_transient(c,qp) );
  libMesh::Real rho_cp = rho*this->_cp(T);
  libMesh::Real T_dot = c.interior_value(this->_T_var, qp);

  libMesh::Real RE_t = rho_cp*T_dot;

  if( this->_enable_thermo_press_calc )
    {
      RE_t -= c.interior_value(this->_p0_var, qp);
    }
  
  return RE_t;
}

template<class Mu, class SH, class TC>
bool GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::element_constraint( bool request_jacobian,
									       libMesh::DiffContext& context,
									       libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  //this->_timer->BeginTimer("LowMachNavierStokesVMSStabilization::element_constraint");
#endif

  //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

#ifdef USE_GRVY_TIMERS
  //this->_timer->EndTimer("LowMachNavierStokesVMSStabilization::element_constraint");
#endif

  return request_jacobian;
}

template<class Mu, class SH, class TC>
bool GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::side_time_derivative( bool request_jacobian,
										 libMesh::DiffContext& context,
										 libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
      //this->_timer->BeginTimer("LowMachNavierStokesVMSStabilization::side_time_derivative");
#endif
      //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

#ifdef USE_GRVY_TIMERS
      //this->_timer->EndTimer("LowMachNavierStokesVMSStabilization::side_time_derivative");
#endif

  return request_jacobian;
}

template<class Mu, class SH, class TC>
bool GRINS::LowMachNavierStokesVMSStabilization<Mu,SH,TC>::side_constraint( bool request_jacobian,
									    libMesh::DiffContext& context,
									    libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  //this->_timer->BeginTimer("LowMachNavierStokesVMSStabilization::side_constraint");
#endif

  //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

#ifdef USE_GRVY_TIMERS
  //this->_timer->EndTimer("LowMachNavierStokesVMSStabilization::side_constraint");
#endif

  return request_jacobian;
}

// Instantiate
template class GRINS::LowMachNavierStokesVMSStabilization<GRINS::ConstantViscosity,GRINS::ConstantSpecificHeat,GRINS::ConstantConductivity>;
