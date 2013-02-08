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
#include "grins/heat_transfer_stab_base.h"

// GRINS
#include "grins/constant_conductivity.h"

namespace GRINS
{
  template<class Conductivity>
  HeatTransferStabilizationBase<Conductivity>::HeatTransferStabilizationBase( const std::string& physics_name, 
									      const GetPot& input )
    : HeatTransferBase<Conductivity>(physics_name,input),
      _stab_helper( input )
  {
    return;
  }

  template<class Conductivity>
  HeatTransferStabilizationBase<Conductivity>::~HeatTransferStabilizationBase()
  {
    return;
  }

  template<class Conductivity>
  void HeatTransferStabilizationBase<Conductivity>::init_context( libMesh::FEMContext& context )
  {
    // First call base class
    HeatTransferBase<Conductivity>::init_context(context);

    // We also need second derivatives, so initialize those.
    context.element_fe_var[this->_T_var]->get_d2phi();

    return;
  }

  template<class Conductivity>
  libMesh::Real HeatTransferStabilizationBase<Conductivity>::compute_res_steady( libMesh::FEMContext& context,
										 unsigned int qp ) const
  {
    libMesh::Gradient grad_T = context.fixed_interior_gradient(this->_T_var, qp);
    libMesh::Tensor hess_T = context.fixed_interior_hessian(this->_T_var, qp);

    libMesh::Real T = context.interior_value(this->_T_var, qp );
    libMesh::Real k = this->_k(T);

    libMesh::RealGradient rhocpU( this->_rho*this->_Cp*context.fixed_interior_value(this->_u_var, qp), 
				  this->_rho*this->_Cp*context.fixed_interior_value(this->_v_var, qp) );
    if(this->_dim == 3)
      rhocpU(2) = this->_rho*this->_Cp*context.fixed_interior_value(this->_w_var, qp);

    return rhocpU*grad_T - k*(hess_T(0,0) + hess_T(1,1) + hess_T(2,2));
  }

  template<class Conductivity>
  libMesh::Real HeatTransferStabilizationBase<Conductivity>::compute_res_transient( libMesh::FEMContext& context,
										    unsigned int qp ) const
  {
    libMesh::Real T_dot = context.interior_value(this->_T_var, qp);

    return this->_rho*this->_Cp*T_dot;
  }

  // Instantiate
  template class HeatTransferStabilizationBase<GRINS::ConstantConductivity>;

} // namespace GRINS
