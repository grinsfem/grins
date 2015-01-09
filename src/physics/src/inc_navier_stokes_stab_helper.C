//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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

//libMesh
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  IncompressibleNavierStokesStabilizationHelper::IncompressibleNavierStokesStabilizationHelper(const GetPot& input)
    : StabilizationHelper(),
      _C( input("Stabilization/tau_constant_vel", input("Stabilization/tau_constant", 1 ) ) ),
      _tau_factor( input("Stabilization/tau_factor_vel", input("Stabilization/tau_factor", 0.5 ) ) ),
      _flow_vars(input)
  {
    return;
  }

  IncompressibleNavierStokesStabilizationHelper::~IncompressibleNavierStokesStabilizationHelper()
  {
    return;
  }

  void IncompressibleNavierStokesStabilizationHelper::init( libMesh::FEMSystem& system )
  {
    _flow_vars.init(&system);
    
    return;
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

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::div_GradU( libMesh::RealTensor& hess_u, 
                                                                                  libMesh::RealTensor& hess_v ) const
  {
    return libMesh::RealGradient( hess_u(0,0) + hess_u(1,1),
                                  hess_v(0,0) + hess_v(1,1) );
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::div_GradU( libMesh::RealTensor& hess_u, 
                                                                                  libMesh::RealTensor& hess_v,
                                                                                  libMesh::RealTensor& hess_w ) const
  {
    return libMesh::RealGradient( hess_u(0,0) + hess_u(1,1) + hess_u(2,2),
                                  hess_v(0,0) + hess_v(1,1) + hess_v(2,2),
                                  hess_w(0,0) + hess_w(1,1) + hess_w(2,2) );
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::div_GradU_T( libMesh::RealTensor& hess_u, 
                                                                                    libMesh::RealTensor& hess_v ) const
  {
    return libMesh::RealGradient( hess_u(0,0) + hess_v(0,1),
                                  hess_u(1,0) + hess_v(1,1) );
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::div_GradU_T( libMesh::RealTensor& hess_u, 
                                                                                    libMesh::RealTensor& hess_v,
                                                                                    libMesh::RealTensor& hess_w ) const
  {
    return libMesh::RealGradient( hess_u(0,0) + hess_v(0,1) + hess_w(0,2),
                                  hess_u(1,0) + hess_v(1,1) + hess_w(1,2),
                                  hess_u(2,0) + hess_v(2,1) + hess_w(2,2) );
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::div_divU_I( libMesh::RealTensor& hess_u, 
                                                                                   libMesh::RealTensor& hess_v ) const
  {
    return libMesh::RealGradient( hess_u(0,0) + hess_v(1,0),
                                  hess_u(0,1) + hess_v(1,1) );
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

    grad_u = context.fixed_interior_gradient(this->_flow_vars.u_var(), qp);
    grad_v = context.fixed_interior_gradient(this->_flow_vars.v_var(), qp);

    libMesh::Real divU = grad_u(0) + grad_v(1);

    if( context.get_system().get_mesh().mesh_dimension() == 3 )
      {
        divU += (context.fixed_interior_gradient(this->_flow_vars.w_var(), qp))(2);
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
    libMesh::RealGradient grad_u, grad_v;

    grad_u = context.fixed_interior_gradient(this->_flow_vars.u_var(), qp);
    grad_v = context.fixed_interior_gradient(this->_flow_vars.v_var(), qp);

    libMesh::Real divU = grad_u(0) + grad_v(1);
    d_res_C_dgradU = 0;
    d_res_C_dgradU(0,0) = 1;
    d_res_C_dgradU(1,1) = 1;

    if( context.get_system().get_mesh().mesh_dimension() == 3 )
      {
        divU += (context.fixed_interior_gradient(this->_flow_vars.w_var(), qp))(2);
        d_res_C_dgradU(2,2) = 1;
      }

    res_C = divU;
  }

  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::compute_res_momentum_steady( AssemblyContext& context,
                                                                                                    unsigned int qp, const libMesh::Real rho, const libMesh::Real mu ) const
  {
    libMesh::RealGradient U( context.fixed_interior_value(this->_flow_vars.u_var(), qp),
                             context.fixed_interior_value(this->_flow_vars.v_var(), qp) );
    if(context.get_system().get_mesh().mesh_dimension() == 3)
      U(2) = context.fixed_interior_value(this->_flow_vars.w_var(), qp);

    libMesh::RealGradient grad_p = context.fixed_interior_gradient(this->_flow_vars.p_var(), qp);

    libMesh::RealGradient grad_u = context.fixed_interior_gradient(this->_flow_vars.u_var(), qp);
    libMesh::RealGradient grad_v = context.fixed_interior_gradient(this->_flow_vars.v_var(), qp);

    libMesh::RealTensor hess_u = context.fixed_interior_hessian(this->_flow_vars.u_var(), qp);
    libMesh::RealTensor hess_v = context.fixed_interior_hessian(this->_flow_vars.v_var(), qp);

    libMesh::RealGradient rhoUdotGradU;
    libMesh::RealGradient divGradU;

    if( context.get_system().get_mesh().mesh_dimension() < 3 )
      {
        rhoUdotGradU = rho*this->UdotGradU( U, grad_u, grad_v );
        divGradU  = this->div_GradU( hess_u, hess_v );
      }
    else
      {
        libMesh::RealGradient grad_w = context.fixed_interior_gradient(this->_flow_vars.w_var(), qp);
        libMesh::RealTensor hess_w = context.fixed_interior_hessian(this->_flow_vars.w_var(), qp);

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
    libMesh::RealGradient U( context.fixed_interior_value(this->_flow_vars.u_var(), qp),
                             context.fixed_interior_value(this->_flow_vars.v_var(), qp) );
    if(context.get_system().get_mesh().mesh_dimension() == 3)
      U(2) = context.fixed_interior_value(this->_flow_vars.w_var(), qp);

    libMesh::RealGradient grad_p = context.fixed_interior_gradient(this->_flow_vars.p_var(), qp);

    libMesh::RealGradient grad_u = context.fixed_interior_gradient(this->_flow_vars.u_var(), qp);
    libMesh::RealGradient grad_v = context.fixed_interior_gradient(this->_flow_vars.v_var(), qp);

    libMesh::RealTensor hess_u = context.fixed_interior_hessian(this->_flow_vars.u_var(), qp);
    libMesh::RealTensor hess_v = context.fixed_interior_hessian(this->_flow_vars.v_var(), qp);

    libMesh::RealGradient rhoUdotGradU;
    libMesh::RealGradient divGradU;

    if( context.get_system().get_mesh().mesh_dimension() < 3 )
      {
        rhoUdotGradU = rho*this->UdotGradU( U, grad_u, grad_v );
        divGradU  = this->div_GradU( hess_u, hess_v );
      }
    else
      {
        libMesh::RealGradient grad_w = context.fixed_interior_gradient(this->_flow_vars.w_var(), qp);
        libMesh::RealTensor hess_w = context.fixed_interior_hessian(this->_flow_vars.w_var(), qp);

        rhoUdotGradU = rho*this->UdotGradU( U, grad_u, grad_v, grad_w );

        divGradU  = this->div_GradU( hess_u, hess_v, hess_w );

        d_res_M_dU(0,2) = rho * grad_u(2);
        d_res_M_dU(1,2) = rho * grad_v(2);

        d_res_Muvw_dgraduvw(2) = rho * U(2);
        d_res_M_dgradp(2,2) = 1;
        d_res_Muvw_dhessuvw(2,2) = mu;

        d_res_M_dU(2,0) = rho * grad_w(0);
        d_res_M_dU(2,1) = rho * grad_w(1);
        d_res_M_dU(2,2) = rho * grad_w(2);
      }

    res_M = rhoUdotGradU + grad_p - mu*divGradU;

    d_res_M_dU(0,0) = rho * grad_u(0);
    d_res_M_dU(0,1) = rho * grad_u(1);
    d_res_M_dU(1,0) = rho * grad_v(0);
    d_res_M_dU(1,1) = rho * grad_v(1);

    d_res_Muvw_dgraduvw(0) = rho * U(0);
    d_res_Muvw_dgraduvw(1) = rho * U(1);

    d_res_M_dgradp(0,0) = 1;
    d_res_M_dgradp(1,1) = 1;

    d_res_Muvw_dhessuvw(0,0) = mu;
    d_res_Muvw_dhessuvw(1,1) = mu;
  }


  libMesh::RealGradient IncompressibleNavierStokesStabilizationHelper::compute_res_momentum_transient( AssemblyContext& context, unsigned int qp, const libMesh::Real rho ) const
  {
    libMesh::RealGradient u_dot( context.interior_rate(this->_flow_vars.u_var(), qp), context.interior_rate(this->_flow_vars.v_var(), qp) );

    if(context.get_system().get_mesh().mesh_dimension() == 3)
      u_dot(2) = context.interior_rate(this->_flow_vars.w_var(), qp);

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
    libMesh::RealGradient u_dot( context.interior_rate(this->_flow_vars.u_var(), qp), context.interior_rate(this->_flow_vars.v_var(), qp) );

    if(context.get_system().get_mesh().mesh_dimension() == 3)
      u_dot(2) = context.interior_rate(this->_flow_vars.w_var(), qp);

    res_M = rho*u_dot;
    d_res_Muvw_duvw = rho;
  }

} // namespace GRINS
