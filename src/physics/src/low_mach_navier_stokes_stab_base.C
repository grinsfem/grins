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

#include "low_mach_navier_stokes_stab_base.h"

template<class Mu, class SH, class TC>
GRINS::LowMachNavierStokesStabilizationBase<Mu,SH,TC>::LowMachNavierStokesStabilizationBase( const std::string& physics_name, 
											   const GetPot& input )
  : GRINS::LowMachNavierStokesBase<Mu,SH,TC>(physics_name,input),
    _stab_helper( input )
{
  this->read_input_options(input);

  return;
}

template<class Mu, class SH, class TC>
GRINS::LowMachNavierStokesStabilizationBase<Mu,SH,TC>::~LowMachNavierStokesStabilizationBase()
{
  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokesStabilizationBase<Mu,SH,TC>::read_input_options( const GetPot& input )
{
  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokesStabilizationBase<Mu,SH,TC>::init_context( libMesh::DiffContext &context )
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
libMesh::Real GRINS::LowMachNavierStokesStabilizationBase<Mu,SH,TC>::compute_res_continuity_steady( libMesh::FEMContext& c,
												   unsigned int qp ) const
{
  libMesh::Real T = c.fixed_interior_value(this->_T_var, qp);
  libMesh::RealGradient grad_T = c.fixed_interior_gradient(this->_T_var, qp);

  libMesh::RealGradient U( c.fixed_interior_value(this->_u_var, qp),
			   c.fixed_interior_value(this->_v_var, qp) );  

  libMesh::RealGradient grad_u, grad_v;

  grad_u = c.fixed_interior_gradient(this->_u_var, qp);
  grad_v = c.fixed_interior_gradient(this->_v_var, qp);

  libMesh::Real divU = grad_u(0) + grad_v(1);

  if( this->_dim == 3 )
    {
      U(2) = c.fixed_interior_value(this->_w_var, qp);
      divU += (c.fixed_interior_gradient(this->_w_var, qp))(2);
    }

  return divU - (U*grad_T)/T;
}

template<class Mu, class SH, class TC>
libMesh::Real GRINS::LowMachNavierStokesStabilizationBase<Mu,SH,TC>::compute_res_continuity_transient( libMesh::FEMContext& c,
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
libMesh::RealGradient GRINS::LowMachNavierStokesStabilizationBase<Mu,SH,TC>::compute_res_momentum_steady( libMesh::FEMContext& c,
													 unsigned int qp ) const
{
  libMesh::Real T = c.fixed_interior_value(this->_T_var, qp);

  libMesh::Real rho = this->compute_rho(T, this->get_p0_transient(c,qp) );

  libMesh::RealGradient U( c.fixed_interior_value(this->_u_var, qp), 
			   c.fixed_interior_value(this->_v_var, qp) );
  if(this->_dim == 3)
    U(2) = c.fixed_interior_value(this->_w_var, qp);

  libMesh::RealGradient grad_p = c.fixed_interior_gradient(this->_p_var, qp);

  libMesh::RealGradient grad_u = c.fixed_interior_gradient(this->_u_var, qp);
  libMesh::RealGradient grad_v = c.fixed_interior_gradient(this->_v_var, qp);

  libMesh::RealTensor hess_u = c.fixed_interior_hessian(this->_u_var, qp);
  libMesh::RealTensor hess_v = c.fixed_interior_hessian(this->_v_var, qp);

  libMesh::RealGradient rhoUdotGradU;
  libMesh::RealGradient divGradU;
  libMesh::RealGradient divGradUT;
  libMesh::RealGradient divdivU;

  if( this->_dim < 3 )
    {
      rhoUdotGradU = rho*_stab_helper.UdotGradU( U, grad_u, grad_v );
      divGradU  = _stab_helper.div_GradU( hess_u, hess_v );
      divGradUT = _stab_helper.div_GradU_T( hess_u, hess_v );
      divdivU   = _stab_helper.div_divU_I( hess_u, hess_v );
    }
  else
    {
      libMesh::RealGradient grad_w = c.fixed_interior_gradient(this->_w_var, qp);
      libMesh::RealTensor hess_w = c.fixed_interior_hessian(this->_w_var, qp);
      
      rhoUdotGradU = rho*_stab_helper.UdotGradU( U, grad_u, grad_v, grad_w );

      divGradU  = _stab_helper.div_GradU( hess_u, hess_v, hess_w );
      divGradUT = _stab_helper.div_GradU_T( hess_u, hess_v, hess_w );
      divdivU   = _stab_helper.div_divU_I( hess_u, hess_v, hess_w );
    }

  libMesh::RealGradient divT = this->_mu(T)*(divGradU + divGradUT - 2.0/3.0*divdivU);

  if( this->_mu.deriv(T) != 0.0 )
    {
      libMesh::Gradient grad_T = c.fixed_interior_gradient(this->_T_var, qp);

      libMesh::Gradient grad_u = c.fixed_interior_gradient(this->_u_var, qp);
      libMesh::Gradient grad_v = c.fixed_interior_gradient(this->_v_var, qp);

      libMesh::Gradient gradTgradu( grad_T*grad_u, grad_T*grad_v );

      libMesh::Gradient gradTgraduT( grad_T(0)*grad_u(0) + grad_T(1)*grad_u(1),
				     grad_T(0)*grad_v(0) + grad_T(1)*grad_v(1) );

      libMesh::Real divU = grad_u(0) + grad_v(1);

      libMesh::Gradient gradTdivU( grad_T(0)*divU, grad_T(1)*divU );

      if(this->_dim == 3)
	{
	  libMesh::Gradient grad_w = c.fixed_interior_gradient(this->_w_var, qp);

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
libMesh::RealGradient GRINS::LowMachNavierStokesStabilizationBase<Mu,SH,TC>::compute_res_momentum_transient( libMesh::FEMContext& c,
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
libMesh::Real GRINS::LowMachNavierStokesStabilizationBase<Mu,SH,TC>::compute_res_energy_steady( libMesh::FEMContext& c,
											       unsigned int qp ) const
{
  libMesh::Real T = c.fixed_interior_value(this->_T_var, qp);
  libMesh::Gradient grad_T = c.fixed_interior_gradient(this->_T_var, qp);
  libMesh::Tensor hess_T = c.fixed_interior_hessian(this->_T_var, qp);

  libMesh::Real rho = this->compute_rho(T, this->get_p0_transient(c,qp) );
  libMesh::Real rho_cp = rho*this->_cp(T);

  libMesh::RealGradient rhocpU( rho_cp*c.fixed_interior_value(this->_u_var, qp), 
				rho_cp*c.fixed_interior_value(this->_v_var, qp) );
  if(this->_dim == 3)
    rhocpU(2) = rho_cp*c.fixed_interior_value(this->_w_var, qp);

  return rhocpU*grad_T - this->_k.deriv(T)*(grad_T*grad_T) - this->_k(T)*(hess_T(0,0) + hess_T(1,1) + hess_T(2,2));
}


template<class Mu, class SH, class TC>
libMesh::Real GRINS::LowMachNavierStokesStabilizationBase<Mu,SH,TC>::compute_res_energy_transient( libMesh::FEMContext& c,
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

// Instantiate
template class GRINS::LowMachNavierStokesStabilizationBase<GRINS::ConstantViscosity,GRINS::ConstantSpecificHeat,GRINS::ConstantConductivity>;
