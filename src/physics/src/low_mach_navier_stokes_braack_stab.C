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

#include "low_mach_navier_stokes_braack_stab.h"

template<class Mu, class SH, class TC>
GRINS::LowMachNavierStokesBraackStabilization<Mu,SH,TC>::LowMachNavierStokesBraackStabilization( const std::string& physics_name, 
											   const GetPot& input )
  : GRINS::LowMachNavierStokesStabilizationBase<Mu,SH,TC>(physics_name,input)
{
  this->read_input_options(input);

  return;
}

template<class Mu, class SH, class TC>
GRINS::LowMachNavierStokesBraackStabilization<Mu,SH,TC>::~LowMachNavierStokesBraackStabilization()
{
  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokesBraackStabilization<Mu,SH,TC>::read_input_options( const GetPot& input )
{
  return;
}

template<class Mu, class SH, class TC>
bool GRINS::LowMachNavierStokesBraackStabilization<Mu,SH,TC>::element_time_derivative( bool request_jacobian,
										    libMesh::DiffContext& context,
										    libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("LowMachNavierStokesBraackStabilization::element_time_derivative");
#endif

  libMesh::FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  this->assemble_continuity_time_deriv( request_jacobian, c, system );
  this->assemble_momentum_time_deriv( request_jacobian, c, system );
  this->assemble_energy_time_deriv( request_jacobian, c, system );

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("LowMachNavierStokesBraackStabilization::element_time_derivative");
#endif
  return request_jacobian;
}

template<class Mu, class SH, class TC>
bool GRINS::LowMachNavierStokesBraackStabilization<Mu,SH,TC>::mass_residual( bool request_jacobian,
									  libMesh::DiffContext& context,
									  libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("LowMachNavierStokesBraackStabilization::mass_residual");
#endif

  libMesh::FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  this->assemble_continuity_mass_residual( request_jacobian, c, system );
  this->assemble_momentum_mass_residual( request_jacobian, c, system );
  this->assemble_energy_mass_residual( request_jacobian, c, system );

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("LowMachNavierStokesBraackStabilization::mass_residual");
#endif
  return request_jacobian;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokesBraackStabilization<Mu,SH,TC>::assemble_continuity_time_deriv( bool request_jacobian,
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
      libMesh::FEBase* fe = c.element_fe_var[this->_u_var];

      libMesh::RealGradient g = this->_stab_helper.compute_g( fe, c, qp );
      libMesh::RealTensor G = this->_stab_helper.compute_G( fe, c, qp );

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
void GRINS::LowMachNavierStokesBraackStabilization<Mu,SH,TC>::assemble_momentum_time_deriv( bool request_jacobian,
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

      libMesh::FEBase* fe = c.element_fe_var[this->_u_var];

      libMesh::RealGradient g = this->_stab_helper.compute_g( fe, c, qp );
      libMesh::RealTensor G = this->_stab_helper.compute_G( fe, c, qp );

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
		     + tau_M*RM_s(0)*rho*U*u_gradphi[i][qp] 
		     + mu*tau_M*RM_s(0)*(u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) 
					 + u_hessphi[i][qp](0,0) + u_hessphi[i][qp](0,1) 
					 - 2.0/3.0*(u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,0)) ) 
		     )*JxW[qp];

	  Fv(i) += ( tau_C*RC_s*u_gradphi[i][qp](1)
		     + tau_M*RM_s(1)*rho*U*u_gradphi[i][qp]
		     + mu*tau_M*RM_s(1)*(u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) 
					 + u_hessphi[i][qp](1,0) + u_hessphi[i][qp](1,1) 
					 - 2.0/3.0*(u_hessphi[i][qp](0,1) + u_hessphi[i][qp](1,1)) ) 
		     )*JxW[qp];

	  if( this->_dim == 3 )
	    {
	      Fu(i) += mu*tau_M*RM_s(0)*(u_hessphi[i][qp](2,2) + u_hessphi[i][qp](0,2) 
					  - 2.0/3.0*u_hessphi[i][qp](2,0))*JxW[qp];

	      Fv(i) += mu*tau_M*RM_s(1)*(u_hessphi[i][qp](2,2) + u_hessphi[i][qp](1,2)
					  - 2.0/3.0*u_hessphi[i][qp](2,1))*JxW[qp];

	      Fw(i) += ( tau_C*RC_s*u_gradphi[i][qp](2)
			 + tau_M*RM_s(2)*rho*U*u_gradphi[i][qp]
			 + mu*tau_M*RM_s(2)*(u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) + u_hessphi[i][qp](2,2)
					     + u_hessphi[i][qp](2,0) + u_hessphi[i][qp](2,1) + u_hessphi[i][qp](2,2)
					     - 2.0/3.0*(u_hessphi[i][qp](0,2) + u_hessphi[i][qp](1,2) 
							+ u_hessphi[i][qp](2,2)) 
					     ) 
			 )*JxW[qp];
	    }
	}

    }
  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokesBraackStabilization<Mu,SH,TC>::assemble_energy_time_deriv( bool request_jacobian,
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

      libMesh::FEBase* fe = c.element_fe_var[this->_u_var];

      libMesh::RealGradient g = this->_stab_helper.compute_g( fe, c, qp );
      libMesh::RealTensor G = this->_stab_helper.compute_G( fe, c, qp );

      libMesh::Real tau_M = this->compute_tau_momentum( c, qp, g, G, rho, U, T, is_steady );
      libMesh::Real tau_E = this->compute_tau_energy( c, qp, g, G, rho, U, T, is_steady );

      libMesh::Real RE_s = this->compute_res_energy_steady( c, qp );
      libMesh::RealGradient RM_s = this->compute_res_momentum_steady( c, qp );

      //std::cout << "tau_E = " << tau_E << ", RE_s = " << RE_s << std::endl;

      libMesh::Real k = this->_k(T);

      for (unsigned int i=0; i != n_T_dofs; i++)
        {
          FT(i) += ( rho_cp*tau_E*RE_s*U*T_gradphi[i][qp]
		     + tau_E*RE_s*k*(T_hessphi[i][qp](0,0) + T_hessphi[i][qp](1,1) + T_hessphi[i][qp](2,2) ) 
		     )*JxW[qp];
	}

    }

  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokesBraackStabilization<Mu,SH,TC>::assemble_continuity_mass_residual( bool request_jacobian,
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
      libMesh::FEBase* fe = c.element_fe_var[this->_u_var];

      libMesh::RealGradient g = this->_stab_helper.compute_g( fe, c, qp );
      libMesh::RealTensor G = this->_stab_helper.compute_G( fe, c, qp );

      libMesh::Real T = c.fixed_interior_value( this->_T_var, qp );
      libMesh::Real rho = this->compute_rho( T, this->get_p0_transient( c, qp ) );

      libMesh::RealGradient U( c.fixed_interior_value( this->_u_var, qp ),
			       c.fixed_interior_value( this->_v_var, qp ) );
      if( this->_dim == 3 )
	U(2) = c.fixed_interior_value( this->_w_var, qp );

      libMesh::Real tau_M = this->compute_tau_momentum( c, qp, g, G, rho, U, T, false );
      libMesh::RealGradient RM_t = this->compute_res_momentum_transient( c, qp );

      libMesh::Real tau_E = this->compute_tau_energy( c, qp, g, G, rho, U, T, false );
      libMesh::Real RE_t = this->compute_res_energy_transient( c, qp );

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) += (tau_E*RE_t*(U*p_dphi[i][qp])/T*JxW[qp]
		    - tau_M*RM_t*p_dphi[i][qp] 
		    )*JxW[qp];
	}
    }

  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokesBraackStabilization<Mu,SH,TC>::assemble_momentum_mass_residual( bool request_jacobian,
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

      libMesh::FEBase* fe = c.element_fe_var[this->_u_var];

      libMesh::RealGradient g = this->_stab_helper.compute_g( fe, c, qp );
      libMesh::RealTensor G = this->_stab_helper.compute_G( fe, c, qp );

      libMesh::Real tau_M = this->compute_tau_momentum( c, qp, g, G, rho, U, T, false );
      libMesh::Real tau_C = this->compute_tau_continuity( tau_M, g, G, U, rho );

      libMesh::Real RC_t = this->compute_res_continuity_transient( c, qp );
      libMesh::RealGradient RM_s = this->compute_res_momentum_steady( c, qp );
      libMesh::RealGradient RM_t = this->compute_res_momentum_transient( c, qp );
      
      libMesh::Real mu = this->_mu(T);

      for (unsigned int i=0; i != n_u_dofs; i++)
        {
	  Fu(i) += ( - tau_C*RC_t*u_gradphi[i][qp](0)
		     - tau_M*RM_t(0)*rho*U*u_gradphi[i][qp]
		     - mu*tau_M*RM_t(0)*(u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) 
					 + u_hessphi[i][qp](0,0) + u_hessphi[i][qp](0,1) 
					 - 2.0/3.0*(u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,0)) ) 
		     )*JxW[qp];

	  Fv(i) += ( -tau_C*RC_t*u_gradphi[i][qp](1)
		     - tau_M*RM_t(1)*rho*U*u_gradphi[i][qp]
		     - mu*tau_M*RM_t(1)*(u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) 
					 + u_hessphi[i][qp](1,0) + u_hessphi[i][qp](1,1) 
					 - 2.0/3.0*(u_hessphi[i][qp](0,1) + u_hessphi[i][qp](1,1)) ) 
		     )*JxW[qp];

	  if( this->_dim == 3 )
	    {
	      Fw(i) -= mu*tau_M*RM_t(0)*(u_hessphi[i][qp](2,2) + u_hessphi[i][qp](0,2) 
					  - 2.0/3.0*u_hessphi[i][qp](2,0))*JxW[qp];

	      Fv(i) -= mu*tau_M*RM_t(1)*(u_hessphi[i][qp](2,2) + u_hessphi[i][qp](1,2)
					  - 2.0/3.0*u_hessphi[i][qp](2,1))*JxW[qp];

	      Fw(i) -= ( tau_C*RC_t*u_gradphi[i][qp](2)
			 + tau_M*RM_s(2)*rho*U*u_gradphi[i][qp]
			 + mu*tau_M*RM_s(2)*(u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) + u_hessphi[i][qp](2,2)
					     + u_hessphi[i][qp](2,0) + u_hessphi[i][qp](2,1) + u_hessphi[i][qp](2,2)
					     - 2.0/3.0*(u_hessphi[i][qp](0,2) + u_hessphi[i][qp](1,2) 
							+ u_hessphi[i][qp](2,2)) 
					     ) 
			 )*JxW[qp];
	    }
	}

    }
  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokesBraackStabilization<Mu,SH,TC>::assemble_energy_mass_residual( bool request_jacobian,
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

      libMesh::FEBase* fe = c.element_fe_var[this->_u_var];

      libMesh::RealGradient g = this->_stab_helper.compute_g( fe, c, qp );
      libMesh::RealTensor G = this->_stab_helper.compute_G( fe, c, qp );

      libMesh::Real tau_M = this->compute_tau_momentum( c, qp, g, G, rho, U, T, false );
      libMesh::Real tau_E = this->compute_tau_energy( c, qp, g, G, rho, U, T, false );

      libMesh::Real RE_s = this->compute_res_energy_steady( c, qp );
      libMesh::Real RE_t = this->compute_res_energy_transient( c, qp );

      libMesh::RealGradient RM_s = this->compute_res_momentum_steady( c, qp );
      libMesh::RealGradient RM_t = this->compute_res_momentum_transient( c, qp );

      libMesh::Real k = this->_k(T);

      for (unsigned int i=0; i != n_T_dofs; i++)
        {
          FT(i) += ( -rho_cp*tau_E*RE_t*U*T_gradphi[i][qp]
		     - tau_E*RE_t*k*(T_hessphi[i][qp](0,0) + T_hessphi[i][qp](1,1) + T_hessphi[i][qp](2,2)) 
		     )*JxW[qp];
	}

    }

  return;
}

template<class Mu, class SH, class TC>
bool GRINS::LowMachNavierStokesBraackStabilization<Mu,SH,TC>::element_constraint( bool request_jacobian,
									       libMesh::DiffContext& context,
									       libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  //this->_timer->BeginTimer("LowMachNavierStokesBraackStabilization::element_constraint");
#endif

  //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

#ifdef USE_GRVY_TIMERS
  //this->_timer->EndTimer("LowMachNavierStokesBraackStabilization::element_constraint");
#endif

  return request_jacobian;
}

template<class Mu, class SH, class TC>
bool GRINS::LowMachNavierStokesBraackStabilization<Mu,SH,TC>::side_time_derivative( bool request_jacobian,
										 libMesh::DiffContext& context,
										 libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
      //this->_timer->BeginTimer("LowMachNavierStokesBraackStabilization::side_time_derivative");
#endif
      //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

#ifdef USE_GRVY_TIMERS
      //this->_timer->EndTimer("LowMachNavierStokesBraackStabilization::side_time_derivative");
#endif

  return request_jacobian;
}

template<class Mu, class SH, class TC>
bool GRINS::LowMachNavierStokesBraackStabilization<Mu,SH,TC>::side_constraint( bool request_jacobian,
									    libMesh::DiffContext& context,
									    libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  //this->_timer->BeginTimer("LowMachNavierStokesBraackStabilization::side_constraint");
#endif

  //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

#ifdef USE_GRVY_TIMERS
  //this->_timer->EndTimer("LowMachNavierStokesBraackStabilization::side_constraint");
#endif

  return request_jacobian;
}

// Instantiate
template class GRINS::LowMachNavierStokesBraackStabilization<GRINS::ConstantViscosity,GRINS::ConstantSpecificHeat,GRINS::ConstantConductivity>;
