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

#include "heat_transfer_stab_base.h"

GRINS::HeatTransferStabilizationBase::HeatTransferStabilizationBase( const std::string& physics_name, 
								     const GetPot& input )
  : GRINS::HeatTransferBase(physics_name,input),
    _stab_helper( input )
{
  this->read_input_options(input);

  return;
}

GRINS::HeatTransferStabilizationBase::~HeatTransferStabilizationBase()
{
  return;
}

void GRINS::HeatTransferStabilizationBase::read_input_options( const GetPot& )
{
  return;
}

void GRINS::HeatTransferStabilizationBase::init_context( libMesh::DiffContext &context )
{
  // First call base class
  GRINS::HeatTransferBase::init_context(context);

  libMesh::FEMContext &c = libmesh_cast_ref<libMesh::FEMContext&>(context);

  // We also need second derivatives, so initialize those.
  c.element_fe_var[this->_T_var]->get_d2phi();

  return;
}

libMesh::Real GRINS::HeatTransferStabilizationBase::compute_res_steady( libMesh::FEMContext& c,
									unsigned int qp ) const
{
  libMesh::Gradient grad_T = c.fixed_interior_gradient(this->_T_var, qp);
  libMesh::Tensor hess_T = c.fixed_interior_hessian(this->_T_var, qp);

  libMesh::RealGradient rhocpU( _rho*_Cp*c.fixed_interior_value(this->_u_var, qp), 
				_rho*_Cp*c.fixed_interior_value(this->_v_var, qp) );
  if(this->_dim == 3)
    rhocpU(2) = _rho*_Cp*c.fixed_interior_value(this->_w_var, qp);

  return rhocpU*grad_T - _k*(hess_T(0,0) + hess_T(1,1) + hess_T(2,2));
}

libMesh::Real GRINS::HeatTransferStabilizationBase::compute_res_transient( libMesh::FEMContext& c,
									   unsigned int qp ) const
{
  libMesh::Real T_dot = c.interior_value(this->_T_var, qp);

  return _rho*_Cp*T_dot;
}
