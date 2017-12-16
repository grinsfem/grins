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
#include "grins/low_mach_navier_stokes_stab_base.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/constant_viscosity.h"
#include "grins/constant_specific_heat.h"
#include "grins/constant_conductivity.h"

namespace GRINS
{

  template<class Mu, class SH, class TC>
  LowMachNavierStokesStabilizationBase<Mu,SH,TC>::LowMachNavierStokesStabilizationBase( const std::string& physics_name,
                                                                                        const GetPot& input )
    : LowMachNavierStokesBase<Mu,SH,TC>(physics_name,PhysicsNaming::low_mach_navier_stokes(),input),
    _stab_helper( physics_name+"StabHelper", input )
  {}

  template<class Mu, class SH, class TC>
  void LowMachNavierStokesStabilizationBase<Mu,SH,TC>::init_context( AssemblyContext& context )
  {
    // First call base class
    LowMachNavierStokesBase<Mu,SH,TC>::init_context(context);

    // We need pressure derivatives
    context.get_element_fe(this->_press_var.p())->get_dphi();

    // We also need second derivatives, so initialize those.
    context.get_element_fe(this->_flow_vars.u())->get_d2phi();
    context.get_element_fe(this->_temp_vars.T())->get_d2phi();
  }

  template<class Mu, class SH, class TC>
  libMesh::Real LowMachNavierStokesStabilizationBase<Mu,SH,TC>::compute_res_continuity_steady( AssemblyContext& context,
                                                                                               unsigned int qp ) const
  {
    libMesh::Real T = context.fixed_interior_value(this->_temp_vars.T(), qp);
    libMesh::RealGradient grad_T = context.fixed_interior_gradient(this->_temp_vars.T(), qp);

    libMesh::RealGradient U( context.fixed_interior_value(this->_flow_vars.u(), qp),
                             context.fixed_interior_value(this->_flow_vars.v(), qp) );

    libMesh::RealGradient grad_u, grad_v;

    grad_u = context.fixed_interior_gradient(this->_flow_vars.u(), qp);
    grad_v = context.fixed_interior_gradient(this->_flow_vars.v(), qp);

    libMesh::Real divU = grad_u(0) + grad_v(1);


    if( this->_flow_vars.dim() == 3 )
      {
        U(2) = context.fixed_interior_value(this->_flow_vars.w(), qp);
        divU += (context.fixed_interior_gradient(this->_flow_vars.w(), qp))(2);
      }

    return divU - (U*grad_T)/T;
  }

  template<class Mu, class SH, class TC>
  libMesh::Real LowMachNavierStokesStabilizationBase<Mu,SH,TC>::compute_res_continuity_transient( AssemblyContext& context,
                                                                                                  unsigned int qp ) const
  {
    libMesh::Real T = context.fixed_interior_value(this->_temp_vars.T(), qp);
    libMesh::Real T_dot = context.interior_value(this->_temp_vars.T(), qp);

    libMesh::Real RC_t = -T_dot/T;

    if( this->_enable_thermo_press_calc )
      {
        libMesh::Real p0 = context.fixed_interior_value(this->_p0_var->p0(), qp);
        libMesh::Real p0_dot = context.interior_value(this->_p0_var->p0(), qp);

        RC_t += p0_dot/p0;
      }

    return RC_t;
  }

  template<class Mu, class SH, class TC>
  libMesh::RealGradient LowMachNavierStokesStabilizationBase<Mu,SH,TC>::compute_res_momentum_steady( AssemblyContext& context,
                                                                                                     unsigned int qp ) const
  {
    libMesh::Real T = context.fixed_interior_value(this->_temp_vars.T(), qp);

    libMesh::Real rho = this->rho(T, this->get_p0_transient(context,qp) );

    libMesh::RealGradient U( context.fixed_interior_value(this->_flow_vars.u(), qp),
                             context.fixed_interior_value(this->_flow_vars.v(), qp) );

    if(this->_flow_vars.dim() == 3)
      U(2) = context.fixed_interior_value(this->_flow_vars.w(), qp);

    libMesh::RealGradient grad_p = context.fixed_interior_gradient(this->_press_var.p(), qp);

    libMesh::RealGradient grad_u = context.fixed_interior_gradient(this->_flow_vars.u(), qp);
    libMesh::RealGradient grad_v = context.fixed_interior_gradient(this->_flow_vars.v(), qp);

    libMesh::RealTensor hess_u = context.fixed_interior_hessian(this->_flow_vars.u(), qp);
    libMesh::RealTensor hess_v = context.fixed_interior_hessian(this->_flow_vars.v(), qp);

    libMesh::RealGradient rhoUdotGradU;
    libMesh::RealGradient divGradU;
    libMesh::RealGradient divGradUT;
    libMesh::RealGradient divdivU;

    if( this->_flow_vars.dim() < 3 )
      {
        rhoUdotGradU = rho*_stab_helper.UdotGradU( U, grad_u, grad_v );
        divGradU  = _stab_helper.div_GradU( hess_u, hess_v );
        divGradUT = _stab_helper.div_GradU_T( hess_u, hess_v );
        divdivU   = _stab_helper.div_divU_I( hess_u, hess_v );
      }
    else
      {
        libMesh::RealGradient grad_w = context.fixed_interior_gradient(this->_flow_vars.w(), qp);
        libMesh::RealTensor hess_w = context.fixed_interior_hessian(this->_flow_vars.w(), qp);

        rhoUdotGradU = rho*_stab_helper.UdotGradU( U, grad_u, grad_v, grad_w );

        divGradU  = _stab_helper.div_GradU( hess_u, hess_v, hess_w );
        divGradUT = _stab_helper.div_GradU_T( hess_u, hess_v, hess_w );
        divdivU   = _stab_helper.div_divU_I( hess_u, hess_v, hess_w );
      }

    libMesh::RealGradient divT = this->_mu(T)*(divGradU + divGradUT - 2.0/3.0*divdivU);

    if( this->_mu.deriv(T) != 0.0 )
      {
        libMesh::Gradient grad_T = context.fixed_interior_gradient(this->_temp_vars.T(), qp);

        libMesh::Gradient grad_u = context.fixed_interior_gradient(this->_flow_vars.u(), qp);
        libMesh::Gradient grad_v = context.fixed_interior_gradient(this->_flow_vars.v(), qp);

        libMesh::Gradient gradTgradu( grad_T*grad_u, grad_T*grad_v );

        libMesh::Gradient gradTgraduT( grad_T(0)*grad_u(0) + grad_T(1)*grad_u(1),
                                       grad_T(0)*grad_v(0) + grad_T(1)*grad_v(1) );

        libMesh::Real divU = grad_u(0) + grad_v(1);

        libMesh::Gradient gradTdivU( grad_T(0)*divU, grad_T(1)*divU );

        if(this->_flow_vars.dim() == 3)
          {
            libMesh::Gradient grad_w = context.fixed_interior_gradient(this->_flow_vars.w(), qp);

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
    if(this->_flow_vars.dim() == 3)
      rhog(2) = rho*this->_g(2);

    return rhoUdotGradU + grad_p - divT - rhog;
  }

  template<class Mu, class SH, class TC>
  libMesh::RealGradient LowMachNavierStokesStabilizationBase<Mu,SH,TC>::compute_res_momentum_transient( AssemblyContext& context,
                                                                                                        unsigned int qp ) const
  {
    libMesh::Real T = context.fixed_interior_value(this->_temp_vars.T(), qp);
    libMesh::Real rho = this->rho(T, this->get_p0_transient(context,qp) );

    libMesh::RealGradient u_dot( context.interior_value(this->_flow_vars.u(), qp), context.interior_value(this->_flow_vars.v(), qp) );

    if(this->_flow_vars.dim() == 3)
      u_dot(2) = context.interior_value(this->_flow_vars.w(), qp);

    return rho*u_dot;
  }

  template<class Mu, class SH, class TC>
  libMesh::Real LowMachNavierStokesStabilizationBase<Mu,SH,TC>::compute_res_energy_steady( AssemblyContext& context,
                                                                                           unsigned int qp ) const
  {
    libMesh::Real T = context.fixed_interior_value(this->_temp_vars.T(), qp);
    libMesh::Gradient grad_T = context.fixed_interior_gradient(this->_temp_vars.T(), qp);
    libMesh::Tensor hess_T = context.fixed_interior_hessian(this->_temp_vars.T(), qp);

    libMesh::Real rho = this->rho(T, this->get_p0_transient(context,qp) );
    libMesh::Real rho_cp = rho*this->_cp(T);

    libMesh::RealGradient rhocpU( rho_cp*context.fixed_interior_value(this->_flow_vars.u(), qp),
                                  rho_cp*context.fixed_interior_value(this->_flow_vars.v(), qp) );

    if(this->_flow_vars.dim() == 3)
      rhocpU(2) = rho_cp*context.fixed_interior_value(this->_flow_vars.w(), qp);

    libMesh::Real hess_term = hess_T(0,0) + hess_T(1,1);
#if LIBMESH_DIM > 2
    hess_term += hess_T(2,2);
#endif

    return rhocpU*grad_T - this->_k.deriv(T)*(grad_T*grad_T) - this->_k(T)*(hess_term);
  }


  template<class Mu, class SH, class TC>
  libMesh::Real LowMachNavierStokesStabilizationBase<Mu,SH,TC>::compute_res_energy_transient( AssemblyContext& context,
                                                                                              unsigned int qp ) const
  {
    libMesh::Real T = context.fixed_interior_value(this->_temp_vars.T(), qp);
    libMesh::Real rho = this->rho(T, this->get_p0_transient(context,qp) );
    libMesh::Real rho_cp = rho*this->_cp(T);
    libMesh::Real T_dot = context.interior_value(this->_temp_vars.T(), qp);

    libMesh::Real RE_t = rho_cp*T_dot;

    if( this->_enable_thermo_press_calc )
      {
        RE_t -= context.interior_value(this->_p0_var->p0(), qp);
      }

    return RE_t;
  }

} // namespace GRINS

// Instantiate
template class GRINS::LowMachNavierStokesStabilizationBase<GRINS::ConstantViscosity,GRINS::ConstantSpecificHeat,GRINS::ConstantConductivity>;
