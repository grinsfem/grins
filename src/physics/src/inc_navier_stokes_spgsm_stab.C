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
#include "grins/inc_navier_stokes_spgsm_stab.h"

// GRINS
#include "grins/assembly_context.h"

//libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{

  IncompressibleNavierStokesSPGSMStabilization::IncompressibleNavierStokesSPGSMStabilization( const std::string& physics_name, 
                                                                                              const GetPot& input )
    : IncompressibleNavierStokesStabilizationBase(physics_name,input)
  {
    this->read_input_options(input);

    return;
  }

  IncompressibleNavierStokesSPGSMStabilization::~IncompressibleNavierStokesSPGSMStabilization()
  {
    return;
  }

  void IncompressibleNavierStokesSPGSMStabilization::element_time_derivative( bool compute_jacobian,
                                                                              AssemblyContext& context,
                                                                              CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("IncompressibleNavierStokesSPGSMStabilization::element_time_derivative");
#endif

    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u_var()).size();

    // Check number of dofs is same for _flow_vars.u_var(), v_var and w_var.
    libmesh_assert (n_u_dofs == context.get_dof_indices(this->_flow_vars.v_var()).size());
    if (this->_dim == 3)
      {
        libmesh_assert (n_u_dofs == context.get_dof_indices(this->_flow_vars.w_var()).size());
      }

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u_var())->get_JxW();

    // The velocity shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_flow_vars.u_var())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u_var()); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v_var()); // R_{v}
    libMesh::DenseSubVector<libMesh::Number> *Fw = NULL;
    if(this->_dim == 3)
      {
        Fw = &context.get_elem_residual(this->_flow_vars.w_var()); // R_{w}
      }

    libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u_var());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::RealGradient U( context.interior_value( this->_flow_vars.u_var(), qp ),
                                 context.interior_value( this->_flow_vars.v_var(), qp ) );
        if( this->_dim == 3 )
          {
            U(2) = context.interior_value( this->_flow_vars.w_var(), qp );
          }
      
        libMesh::Real tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, this->_rho, U, this->_mu, this->_is_steady );
        libMesh::Real tau_C = this->_stab_helper.compute_tau_continuity( tau_M, g );

        libMesh::RealGradient RM_s = this->_stab_helper.compute_res_momentum_steady( context, qp, _rho, _mu );
        libMesh::Real RC = this->_stab_helper.compute_res_continuity( context, qp );

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) += ( - tau_C*RC*u_gradphi[i][qp](0)
                       - tau_M*RM_s(0)*_rho*U*u_gradphi[i][qp]  )*JxW[qp];

            Fv(i) += ( - tau_C*RC*u_gradphi[i][qp](1)
                       - tau_M*RM_s(1)*_rho*U*u_gradphi[i][qp] )*JxW[qp];

            if( this->_dim == 3 )
              {
                (*Fw)(i) += ( - tau_C*RC*u_gradphi[i][qp](2)
                              - tau_M*RM_s(2)*_rho*U*u_gradphi[i][qp] )*JxW[qp];
              }
          }

        if( compute_jacobian )
          {
            libmesh_not_implemented();
          }

      }

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("IncompressibleNavierStokesSPGSMStabilization::element_time_derivative");
#endif

    return;
  }

  void IncompressibleNavierStokesSPGSMStabilization::element_constraint( bool compute_jacobian,
                                                                         AssemblyContext& context,
                                                                         CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("IncompressibleNavierStokesSPGSMStabilization::element_constraint");
#endif

    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.get_dof_indices(this->_flow_vars.p_var()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u_var())->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
      context.get_element_fe(this->_flow_vars.p_var())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(this->_flow_vars.p_var()); // R_{p}

    libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u_var());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );
      
        libMesh::RealGradient U( context.interior_value( this->_flow_vars.u_var(), qp ),
                                 context.interior_value( this->_flow_vars.v_var(), qp ) );
        if( this->_dim == 3 )
          {
            U(2) = context.interior_value( this->_flow_vars.w_var(), qp );
          }

        libMesh::Real tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, this->_rho, U, this->_mu, this->_is_steady );

        libMesh::RealGradient RM_s = this->_stab_helper.compute_res_momentum_steady( context, qp, _rho, _mu );

        for (unsigned int i=0; i != n_p_dofs; i++)
          {
            Fp(i) += tau_M*RM_s*p_dphi[i][qp]*JxW[qp];
          }

        if( compute_jacobian )
          {
            libmesh_not_implemented();
          }

      }

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("IncompressibleNavierStokesSPGSMStabilization::element_constraint");
#endif

    return;
  }

  void IncompressibleNavierStokesSPGSMStabilization::mass_residual( bool compute_jacobian,
                                                                    AssemblyContext& context,
                                                                    CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("IncompressibleNavierStokesSPGSMStabilization::mass_residual");
#endif

    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.get_dof_indices(this->_flow_vars.p_var()).size();
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u_var()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u_var())->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
      context.get_element_fe(this->_flow_vars.p_var())->get_dphi();

    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_flow_vars.u_var())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u_var()); // R_{p}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v_var()); // R_{p}
    libMesh::DenseSubVector<libMesh::Number> *Fw = NULL;
    if(this->_dim == 3)
      Fw = &context.get_elem_residual(this->_flow_vars.w_var()); // R_{w}

    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(this->_flow_vars.p_var()); // R_{p}

    libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u_var());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::RealGradient U( context.fixed_interior_value( this->_flow_vars.u_var(), qp ),
                                 context.fixed_interior_value( this->_flow_vars.v_var(), qp ) );
        if( this->_dim == 3 )
          {
            U(2) = context.fixed_interior_value( this->_flow_vars.w_var(), qp );
          }
      
        libMesh::Real tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, this->_rho, U, this->_mu, false );

        libMesh::RealGradient RM_t = this->_stab_helper.compute_res_momentum_transient( context, qp, _rho );

        for (unsigned int i=0; i != n_p_dofs; i++)
          {
            Fp(i) += tau_M*RM_t*p_dphi[i][qp]*JxW[qp];
          }

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) += tau_M*RM_t(0)*_rho*U*u_gradphi[i][qp]*JxW[qp];

            Fv(i) += tau_M*RM_t(1)*_rho*U*u_gradphi[i][qp]*JxW[qp];

            if( this->_dim == 3 )
              {
                (*Fw)(i) += tau_M*RM_t(2)*_rho*U*u_gradphi[i][qp]*JxW[qp];
              }
          }

        if( compute_jacobian )
          {
            libmesh_not_implemented();
          }

      }

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("IncompressibleNavierStokesSPGSMStabilization::mass_residual");
#endif

    return;
  }

} // end namespace GRINS
