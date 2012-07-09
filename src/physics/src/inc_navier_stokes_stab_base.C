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

#include "inc_navier_stokes_stab_base.h"

GRINS::IncompressibleNavierStokesStabilizationBase::IncompressibleNavierStokesStabilizationBase( const std::string& physics_name, 
												 const GetPot& input )
  : GRINS::IncompressibleNavierStokesBase(physics_name,input),
    _stab_helper( input )
{
  this->read_input_options(input);

  return;
}

GRINS::IncompressibleNavierStokesStabilizationBase::~IncompressibleNavierStokesStabilizationBase()
{
  return;
}

void GRINS::IncompressibleNavierStokesStabilizationBase::read_input_options( const GetPot& input )
{
  return;
}

void GRINS::IncompressibleNavierStokesStabilizationBase::init_context( libMesh::DiffContext &context )
{
  // First call base class
  GRINS::IncompressibleNavierStokesBase::init_context(context);

  libMesh::FEMContext &c = libmesh_cast_ref<libMesh::FEMContext&>(context);
  
  // We need pressure derivatives
  c.element_fe_var[this->_p_var]->get_dphi();

  // We also need second derivatives, so initialize those.
  c.element_fe_var[this->_u_var]->get_d2phi();

  return;
}

libMesh::Real GRINS::IncompressibleNavierStokesStabilizationBase::compute_res_continuity( libMesh::FEMContext& c,
											  unsigned int qp ) const
{
  libMesh::RealGradient grad_u, grad_v;

  grad_u = c.fixed_interior_gradient(this->_u_var, qp);
  grad_v = c.fixed_interior_gradient(this->_v_var, qp);

  libMesh::Real divU = grad_u(0) + grad_v(1);

  if( this->_dim == 3 )
    {
      divU += (c.fixed_interior_gradient(this->_w_var, qp))(2);
    }

  return divU;
}

libMesh::RealGradient GRINS::IncompressibleNavierStokesStabilizationBase::compute_res_momentum_steady( libMesh::FEMContext& c,
												       unsigned int qp ) const
{
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

  if( this->_dim < 3 )
    {
      rhoUdotGradU = this->_rho*_stab_helper.UdotGradU( U, grad_u, grad_v );
      divGradU  = _stab_helper.div_GradU( hess_u, hess_v );
    }
  else
    {
      libMesh::RealGradient grad_w = c.fixed_interior_gradient(this->_w_var, qp);
      libMesh::RealTensor hess_w = c.fixed_interior_hessian(this->_w_var, qp);
      
      rhoUdotGradU = this->_rho*_stab_helper.UdotGradU( U, grad_u, grad_v, grad_w );

      divGradU  = _stab_helper.div_GradU( hess_u, hess_v, hess_w );
    }

  return rhoUdotGradU + grad_p - this->_mu*divGradU;
}

libMesh::RealGradient GRINS::IncompressibleNavierStokesStabilizationBase::compute_res_momentum_transient( libMesh::FEMContext& c,
													  unsigned int qp ) const
{
  libMesh::RealGradient u_dot( c.interior_value(this->_u_var, qp), c.interior_value(this->_v_var, qp) );

  if(this->_dim == 3)
    u_dot(2) = c.interior_value(this->_w_var, qp);

  return this->_rho*u_dot;
}
