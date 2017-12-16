//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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

// GRINS
#include "grins/variables_parsing.h"
#include "grins/variable_warehouse.h"
#include "grins/multi_component_vector_variable.h"
#include "grins/single_variable.h"
#include "grins/physics_naming.h"

//libMesh
#include "libmesh/getpot.h"
#include "libmesh/system.h"
#include "libmesh/mesh.h"

namespace GRINS
{

  HeatTransferStabilizationHelper::HeatTransferStabilizationHelper
  (const std::string & helper_name,
   const GetPot& input)
    : StabilizationHelper(helper_name),
      _C(1),
      _tau_factor(0.5),
      _temp_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<PrimitiveTempFEVariables>(VariablesParsing::temp_variable_name(input,PhysicsNaming::heat_transfer(),VariablesParsing::PHYSICS))),
      _flow_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<VelocityVariable>(VariablesParsing::velocity_variable_name(input,PhysicsNaming::heat_transfer(),VariablesParsing::PHYSICS))),
      _press_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PressureFEVariable>(VariablesParsing::press_variable_name(input,PhysicsNaming::heat_transfer(),VariablesParsing::PHYSICS)))
  {
    if (input.have_variable("Stabilization/tau_constant_T"))
      this->set_parameter
        (_C, input, "Stabilization/tau_constant_T", _C );
    else
      this->set_parameter
        (_C, input, "Stabilization/tau_constant", _C );

    if (input.have_variable("Stabilization/tau_factor_T"))
      this->set_parameter
        (_tau_factor, input, "Stabilization/tau_factor_T", _tau_factor );
    else
      this->set_parameter
        (_tau_factor, input, "Stabilization/tau_factor", _tau_factor );
  }

  HeatTransferStabilizationHelper::~HeatTransferStabilizationHelper()
  {
    return;
  }

  libMesh::Real HeatTransferStabilizationHelper::compute_res_energy_steady( AssemblyContext& context,
                                                                            unsigned int qp,
                                                                            const libMesh::Real rho,
                                                                            const libMesh::Real Cp,
                                                                            const libMesh::Real k ) const
  {
    libMesh::Gradient grad_T = context.fixed_interior_gradient(this->_temp_vars.T(), qp);
    libMesh::Tensor hess_T = context.fixed_interior_hessian(this->_temp_vars.T(), qp);

    libMesh::RealGradient rhocpU( rho*Cp*context.fixed_interior_value(this->_flow_vars.u(), qp),
                                  rho*Cp*context.fixed_interior_value(this->_flow_vars.v(), qp) );
    if(this->_flow_vars.dim() == 3)
      rhocpU(2) = rho*Cp*context.fixed_interior_value(this->_flow_vars.w(), qp);

    return rhocpU*grad_T - k*(hess_T(0,0) + hess_T(1,1) + hess_T(2,2));
  }

  void HeatTransferStabilizationHelper::compute_res_energy_steady_and_derivs
  ( AssemblyContext& context,
    unsigned int qp,
    const libMesh::Real rho,
    const libMesh::Real Cp,
    const libMesh::Real k,
    libMesh::Real &res,
    libMesh::Real &d_res_dT,
    libMesh::Gradient &d_res_dgradT,
    libMesh::Tensor   &d_res_dhessT,
    libMesh::Gradient &d_res_dU
    ) const
  {
    libMesh::Gradient grad_T = context.fixed_interior_gradient(this->_temp_vars.T(), qp);
    libMesh::Tensor hess_T = context.fixed_interior_hessian(this->_temp_vars.T(), qp);

    libMesh::RealGradient rhocpU( rho*Cp*context.fixed_interior_value(this->_flow_vars.u(), qp),
                                  rho*Cp*context.fixed_interior_value(this->_flow_vars.v(), qp) );
    if(this->_flow_vars.dim() == 3)
      rhocpU(2) = rho*Cp*context.fixed_interior_value(this->_flow_vars.w(), qp);

    res = rhocpU*grad_T - k*(hess_T(0,0) + hess_T(1,1) + hess_T(2,2));
    d_res_dT = 0;
    d_res_dgradT = rhocpU;
    d_res_dhessT = 0;
    d_res_dhessT(0,0) = -k;
    d_res_dhessT(1,1) = -k;
    d_res_dhessT(2,2) = -k;
    d_res_dU = rho * Cp * grad_T;
  }


  libMesh::Real HeatTransferStabilizationHelper::compute_res_energy_transient( AssemblyContext& context,
                                                                               unsigned int qp,
                                                                               const libMesh::Real rho,
                                                                               const libMesh::Real Cp ) const
  {
    libMesh::Real T_dot;
    context.interior_rate(this->_temp_vars.T(), qp, T_dot);

    return rho*Cp*T_dot;
  }


  void HeatTransferStabilizationHelper::compute_res_energy_transient_and_derivs
  ( AssemblyContext& context,
    unsigned int qp,
    const libMesh::Real rho,
    const libMesh::Real Cp,
    libMesh::Real &res,
    libMesh::Real &d_res_dTdot
    ) const
  {
    libMesh::Real T_dot;
    context.interior_rate(this->_temp_vars.T(), qp, T_dot);

    res = rho*Cp*T_dot;
    d_res_dTdot = rho*Cp;
  }

} // namespace GRINS
