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
#include "grins/inc_navier_stokes_stab_helper.h"

// GRINS
#include "grins/variables_parsing.h"
#include "grins/variable_warehouse.h"
#include "grins/multi_component_vector_variable.h"
#include "grins/single_variable.h"
#include "grins/physics_naming.h"

//libMesh
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  IncompressibleNavierStokesStabilizationHelper::IncompressibleNavierStokesStabilizationHelper
  (const std::string & helper_name,
   const GetPot& input)
    : StabilizationHelper(helper_name),
      _C(1),
      _tau_factor(0.5),
      _flow_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<VelocityVariable>(VariablesParsing::velocity_variable_name(input,PhysicsNaming::incompressible_navier_stokes(),VariablesParsing::PHYSICS))),
      _press_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PressureFEVariable>(VariablesParsing::press_variable_name(input,PhysicsNaming::incompressible_navier_stokes(),VariablesParsing::PHYSICS)))
  {
    if (input.have_variable("Stabilization/tau_constant_vel"))
      this->set_parameter
        (_C, input, "Stabilization/tau_constant_vel", _C );
    else
      this->set_parameter
        (_C, input, "Stabilization/tau_constant", _C );

    if (input.have_variable("Stabilization/tau_factor_vel"))
      this->set_parameter
        (_tau_factor, input, "Stabilization/tau_factor_vel", _tau_factor );
    else
      this->set_parameter
        (_tau_factor, input, "Stabilization/tau_factor", _tau_factor );

    return;
  }

  IncompressibleNavierStokesStabilizationHelper::~IncompressibleNavierStokesStabilizationHelper()
  {
    return;
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::UdotGradU( libMesh::Gradient& U,
                                                                                  libMesh::Gradient& grad_u ) const
  {
    return libMesh::RealGradient( U*grad_u );
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::UdotGradU( libMesh::Gradient& U,
                                                                                  libMesh::Gradient& grad_u,
                                                                                  libMesh::Gradient& grad_v ) const
  {
    return libMesh::RealGradient( U*grad_u, U*grad_v );
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::UdotGradU( libMesh::Gradient& U,
                                                                                  libMesh::Gradient& grad_u,
                                                                                  libMesh::Gradient& grad_v,
                                                                                  libMesh::Gradient& grad_w ) const
  {
    return libMesh::RealGradient( U*grad_u, U*grad_v, U*grad_w );
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::div_GradU( libMesh::RealTensor& hess_u ) const
  {
    return libMesh::RealGradient( hess_u(0,0) );
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::div_GradU( libMesh::RealTensor& hess_u,
                                                                                  libMesh::RealTensor& hess_v ) const
  {
    return libMesh::RealGradient( hess_u(0,0) + hess_u(1,1),
                                  hess_v(0,0) + hess_v(1,1) );
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::div_GradU_axi( libMesh::Real r,
                                                                                      const libMesh::Gradient& U,
                                                                                      const libMesh::Gradient& grad_u,
                                                                                      const libMesh::Gradient& grad_v,
                                                                                      const libMesh::RealTensor& hess_u,
                                                                                      const libMesh::RealTensor& hess_v ) const
  {
    return libMesh::RealGradient( hess_u(0,0) + hess_u(1,1) + (grad_u(0) - U(0)/r)/r,
                                  hess_v(0,0) + hess_v(1,1) + grad_v(0)/r );
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::div_GradU( libMesh::RealTensor& hess_u,
                                                                                  libMesh::RealTensor& hess_v,
                                                                                  libMesh::RealTensor& hess_w ) const
  {
    return libMesh::RealGradient( hess_u(0,0) + hess_u(1,1) + hess_u(2,2),
                                  hess_v(0,0) + hess_v(1,1) + hess_v(2,2),
                                  hess_w(0,0) + hess_w(1,1) + hess_w(2,2) );
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::div_GradU_T( libMesh::RealTensor& hess_u ) const
  {
    return libMesh::RealGradient( hess_u(0,0) );
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::div_GradU_T( libMesh::RealTensor& hess_u,
                                                                                    libMesh::RealTensor& hess_v ) const
  {
    return libMesh::RealGradient( hess_u(0,0) + hess_v(0,1),
                                  hess_u(1,0) + hess_v(1,1));
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::div_GradU_T_axi( libMesh::Real r,
                                                                                        const libMesh::Gradient& U,
                                                                                        const libMesh::Gradient& grad_u,
                                                                                        const libMesh::RealTensor& hess_u,
                                                                                        const libMesh::RealTensor& hess_v ) const
  {
    return libMesh::RealGradient( hess_u(0,0) + hess_v(0,1) + (grad_u(0) - U(0)/r)/r,
                                  hess_u(1,0) + hess_v(1,1) + grad_u(1)/r );
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::div_GradU_T( libMesh::RealTensor& hess_u,
                                                                                    libMesh::RealTensor& hess_v,
                                                                                    libMesh::RealTensor& hess_w ) const
  {
    return libMesh::RealGradient( hess_u(0,0) + hess_v(0,1) + hess_w(0,2),
                                  hess_u(1,0) + hess_v(1,1) + hess_w(1,2),
                                  hess_u(2,0) + hess_v(2,1) + hess_w(2,2) );
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::div_divU_I( libMesh::RealTensor& hess_u ) const
  {
    return libMesh::RealGradient( hess_u(0,0) );
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::div_divU_I( libMesh::RealTensor& hess_u,
                                                                                   libMesh::RealTensor& hess_v ) const
  {
    return libMesh::RealGradient( hess_u(0,0) + hess_v(1,0),
                                  hess_u(0,1) + hess_v(1,1) );
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::div_divU_I_axi( libMesh::Real r,
                                                                                       const libMesh::Gradient& U,
                                                                                       const libMesh::Gradient& grad_u,
                                                                                       const libMesh::RealTensor& hess_u,
                                                                                       const libMesh::RealTensor& hess_v ) const
  {
    return libMesh::RealGradient( hess_u(0,0) + hess_v(1,0) + (grad_u(0) - U(0)/r)/r,
                                  hess_u(0,1) + hess_v(1,1) + grad_u(1)/r );
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::div_divU_I( libMesh::RealTensor& hess_u,
                                                                                   libMesh::RealTensor& hess_v,
                                                                                   libMesh::RealTensor& hess_w) const
  {
    return libMesh::RealGradient( hess_u(0,0) + hess_v(1,0) + hess_w(2,0),
                                  hess_u(0,1) + hess_v(1,1) + hess_w(2,1),
                                  hess_u(0,2) + hess_v(1,2) + hess_w(2,2) );
  }

  libMesh::Real IncompressibleNavierStokesStabilizationHelper::compute_res_continuity( AssemblyContext& context,
                                                                                       unsigned int qp ) const
  {
    libMesh::RealGradient grad_u, grad_v;

    grad_u = context.fixed_interior_gradient(this->_flow_vars.u(), qp);

    libMesh::Real divU = grad_u(0);

    if( this->_flow_vars.dim() > 1 )
      {
        divU += (context.fixed_interior_gradient(this->_flow_vars.v(), qp))(1);
      }

    if( this->_flow_vars.dim() == 3 )
      {
        divU += (context.fixed_interior_gradient(this->_flow_vars.w(), qp))(2);
      }

    return divU;
  }

  void IncompressibleNavierStokesStabilizationHelper::compute_res_continuity_and_derivs
  ( AssemblyContext& context,
    unsigned int qp,
    libMesh::Real   &res_C,
    libMesh::Tensor &d_res_C_dgradU
    ) const
  {
    libMesh::RealGradient grad_u =
      context.fixed_interior_gradient(this->_flow_vars.u(), qp);

    libMesh::Real divU = grad_u(0);
    d_res_C_dgradU = 0;
    d_res_C_dgradU(0,0) = 1;

    if( this->_flow_vars.dim() > 1 )
      {
        divU += (context.fixed_interior_gradient(this->_flow_vars.v(), qp))(1);
        d_res_C_dgradU(1,1) = 1;
      }

    if( this->_flow_vars.dim() == 3 )
      {
        divU += (context.fixed_interior_gradient(this->_flow_vars.w(), qp))(2);
        d_res_C_dgradU(2,2) = 1;
      }

    res_C = divU;
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::compute_res_momentum_steady( AssemblyContext& context,
                                                                                                    unsigned int qp, const libMesh::Real rho, const libMesh::Real mu ) const
  {
    libMesh::RealGradient U( context.fixed_interior_value(this->_flow_vars.u(), qp) );
    if(this->_flow_vars.dim() > 1)
      U(1) = context.fixed_interior_value(this->_flow_vars.v(), qp);
    if(this->_flow_vars.dim() == 3)
      U(2) = context.fixed_interior_value(this->_flow_vars.w(), qp);

    libMesh::RealGradient grad_p = context.fixed_interior_gradient(this->_press_var.p(), qp);

    libMesh::RealGradient grad_u = context.fixed_interior_gradient(this->_flow_vars.u(), qp);

    libMesh::RealTensor hess_u = context.fixed_interior_hessian(this->_flow_vars.u(), qp);

    libMesh::RealGradient rhoUdotGradU;
    libMesh::RealGradient divGradU;

    if( this->_flow_vars.dim() == 1 )
      {
        rhoUdotGradU = rho*this->UdotGradU( U, grad_u );
        divGradU  = this->div_GradU( hess_u );
      }
    else if( this->_flow_vars.dim() == 2 )
      {
        libMesh::RealGradient grad_v = context.fixed_interior_gradient(this->_flow_vars.v(), qp);
        libMesh::RealTensor hess_v = context.fixed_interior_hessian(this->_flow_vars.v(), qp);

        rhoUdotGradU = rho*this->UdotGradU( U, grad_u, grad_v );
        divGradU  = this->div_GradU( hess_u, hess_v );
      }
    else
      {
        libMesh::RealGradient grad_v = context.fixed_interior_gradient(this->_flow_vars.v(), qp);
        libMesh::RealTensor hess_v = context.fixed_interior_hessian(this->_flow_vars.v(), qp);

        libMesh::RealGradient grad_w = context.fixed_interior_gradient(this->_flow_vars.w(), qp);
        libMesh::RealTensor hess_w = context.fixed_interior_hessian(this->_flow_vars.w(), qp);

        rhoUdotGradU = rho*this->UdotGradU( U, grad_u, grad_v, grad_w );

        divGradU  = this->div_GradU( hess_u, hess_v, hess_w );
      }

    return rhoUdotGradU + grad_p - mu*divGradU;
  }

  void IncompressibleNavierStokesStabilizationHelper::compute_res_momentum_steady_and_derivs
  ( AssemblyContext& context,
    unsigned int qp, const libMesh::Real rho, const libMesh::Real mu,
    libMesh::Gradient &res_M,
    libMesh::Tensor   &d_res_M_dgradp,
    libMesh::Tensor   &d_res_M_dU,
    libMesh::Gradient &d_res_Muvw_dgraduvw,
    libMesh::Tensor   &d_res_Muvw_dhessuvw
    ) const
  {
    libMesh::RealGradient U( context.fixed_interior_value(this->_flow_vars.u(), qp) );
    if(this->_flow_vars.dim() > 1)
      U(1) = context.fixed_interior_value(this->_flow_vars.v(), qp);
    if(this->_flow_vars.dim() == 3)
      U(2) = context.fixed_interior_value(this->_flow_vars.w(), qp);

    libMesh::RealGradient grad_p = context.fixed_interior_gradient(this->_press_var.p(), qp);

    libMesh::RealGradient grad_u = context.fixed_interior_gradient(this->_flow_vars.u(), qp);

    libMesh::RealTensor hess_u = context.fixed_interior_hessian(this->_flow_vars.u(), qp);

    libMesh::RealGradient rhoUdotGradU;
    libMesh::RealGradient divGradU;

    if( this->_flow_vars.dim() == 1 )
      {
        rhoUdotGradU = rho*this->UdotGradU( U, grad_u );
        divGradU  = this->div_GradU( hess_u );
      }
    else if( this->_flow_vars.dim() == 2 )
      {
        libMesh::RealGradient grad_v = context.fixed_interior_gradient(this->_flow_vars.v(), qp);
        libMesh::RealTensor hess_v = context.fixed_interior_hessian(this->_flow_vars.v(), qp);

        d_res_M_dU(0,1) = rho * grad_u(1);
        d_res_M_dU(1,0) = rho * grad_v(0);
        d_res_M_dU(1,1) = rho * grad_v(1);

        d_res_Muvw_dgraduvw(1) = rho * U(1);
        d_res_M_dgradp(1,1) = 1;
        d_res_Muvw_dhessuvw(1,1) = mu;

        if( this->_flow_vars.dim() == 3 )
          {
            libMesh::RealGradient grad_w = context.fixed_interior_gradient(this->_flow_vars.w(), qp);
            libMesh::RealTensor hess_w = context.fixed_interior_hessian(this->_flow_vars.w(), qp);

            rhoUdotGradU = rho*this->UdotGradU( U, grad_u, grad_v, grad_w );

            divGradU  = this->div_GradU( hess_u, hess_v, hess_w );

            d_res_M_dU(0,2) = rho * grad_u(2);
            d_res_M_dU(1,2) = rho * grad_v(2);
            d_res_M_dU(2,0) = rho * grad_w(0);
            d_res_M_dU(2,1) = rho * grad_w(1);
            d_res_M_dU(2,2) = rho * grad_w(2);

            d_res_Muvw_dgraduvw(2) = rho * U(2);
            d_res_M_dgradp(2,2) = 1;
            d_res_Muvw_dhessuvw(2,2) = mu;
          }
        else
          {
            rhoUdotGradU = rho*this->UdotGradU( U, grad_u, grad_v );
            divGradU  = this->div_GradU( hess_u, hess_v );
          }
      }

    res_M = rhoUdotGradU + grad_p - mu*divGradU;

    d_res_M_dU(0,0) = rho * grad_u(0);

    d_res_Muvw_dgraduvw(0) = rho * U(0);

    d_res_M_dgradp(0,0) = 1;

    d_res_Muvw_dhessuvw(0,0) = mu;
  }


  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::compute_res_momentum_transient( AssemblyContext& context, unsigned int qp, const libMesh::Real rho ) const
  {
    libMesh::RealGradient u_dot;
    context.interior_rate(this->_flow_vars.u(), qp, u_dot(0));
    if(this->_flow_vars.dim() > 1)
      context.interior_rate(this->_flow_vars.v(), qp, u_dot(1));
    if(this->_flow_vars.dim() == 3)
      context.interior_rate(this->_flow_vars.w(), qp, u_dot(2));

    return rho*u_dot;
  }


  void IncompressibleNavierStokesStabilizationHelper::compute_res_momentum_transient_and_derivs
  ( AssemblyContext& context,
    unsigned int qp,
    const libMesh::Real rho,
    libMesh::RealGradient &res_M,
    libMesh::Real &d_res_Muvw_duvw
    ) const
  {
    libMesh::RealGradient u_dot;
    context.interior_rate(this->_flow_vars.u(), qp, u_dot(0));
    if(this->_flow_vars.dim() > 1)
      context.interior_rate(this->_flow_vars.v(), qp, u_dot(1));
    if(this->_flow_vars.dim() == 3)
      context.interior_rate(this->_flow_vars.w(), qp, u_dot(2));

    res_M = rho*u_dot;
    d_res_Muvw_duvw = rho;
  }

} // namespace GRINS
