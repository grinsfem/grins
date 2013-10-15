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


// This class
#include "grins/heat_transfer_stab_helper.h"

//libMesh
#include "libmesh/getpot.h"
#include "libmesh/system.h"
#include "libmesh/mesh.h"

namespace GRINS
{

  HeatTransferStabilizationHelper::HeatTransferStabilizationHelper(const GetPot& input)
    : StabilizationHelper(),
      _C( input("Stabilization/tau_constant_T", input("Stabilization/tau_constant", 1 ) ) ),
      _tau_factor( input("Stabilization/tau_factor_T", input("Stabilization/tau_factor", 0.5 ) ) ),
      _temp_vars(input),
      _flow_vars(input)
  {
    return;
  }

  HeatTransferStabilizationHelper::~HeatTransferStabilizationHelper()
  {
    return;
  }

  void HeatTransferStabilizationHelper::init( libMesh::FEMSystem& system )
  {
    _temp_vars.init(&system);
    _flow_vars.init(&system);

    return;
  }

  libMesh::Real HeatTransferStabilizationHelper::compute_res_energy_steady( AssemblyContext& context,
                                                                            unsigned int qp,
                                                                            const libMesh::Real rho,
                                                                            const libMesh::Real Cp,
                                                                            const libMesh::Real k ) const
  {
    libMesh::Gradient grad_T = context.fixed_interior_gradient(this->_temp_vars.T_var(), qp);
    libMesh::Tensor hess_T = context.fixed_interior_hessian(this->_temp_vars.T_var(), qp);

    libMesh::RealGradient rhocpU( rho*Cp*context.fixed_interior_value(this->_flow_vars.u_var(), qp), 
                                  rho*Cp*context.fixed_interior_value(this->_flow_vars.v_var(), qp) );
    if(context.get_system().get_mesh().mesh_dimension() == 3)
      rhocpU(2) = rho*Cp*context.fixed_interior_value(this->_flow_vars.w_var(), qp);

    return rhocpU*grad_T - k*(hess_T(0,0) + hess_T(1,1) + hess_T(2,2));
  }

  libMesh::Real HeatTransferStabilizationHelper::compute_res_energy_transient( AssemblyContext& context,
                                                                               unsigned int qp,
                                                                               const libMesh::Real rho,
                                                                               const libMesh::Real Cp ) const
  {
    libMesh::Real T_dot = context.interior_value(this->_temp_vars.T_var(), qp);

    return rho*Cp*T_dot;
  }

} // namespace GRINS
