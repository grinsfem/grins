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
#include "grins/inc_navier_stokes_spgsm_stab.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/inc_nav_stokes_macro.h"

//libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{

  template<class Mu>
  IncompressibleNavierStokesSPGSMStabilization<Mu>::IncompressibleNavierStokesSPGSMStabilization( const std::string& physics_name,
                                                                                                  const GetPot& input )
    : IncompressibleNavierStokesStabilizationBase<Mu>(physics_name,input)
  {}

  template<class Mu>
  void IncompressibleNavierStokesSPGSMStabilization<Mu>::element_time_derivative
  ( bool compute_jacobian,
    AssemblyContext & context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();

    // Check number of dofs is same for _flow_vars.u(), v_var and w_var.
    libmesh_assert (n_u_dofs == context.get_dof_indices(this->_flow_vars.v()).size());


    if (this->_flow_vars.dim() == 3)
      libmesh_assert (n_u_dofs == context.get_dof_indices(this->_flow_vars.w()).size());

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // The velocity shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_flow_vars.u())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u()); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v()); // R_{v}
    libMesh::DenseSubVector<libMesh::Number> *Fw = NULL;
    if(this->_flow_vars.dim() == 3)
      Fw = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}

    libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::RealGradient U( context.interior_value( this->_flow_vars.u(), qp ),
                                 context.interior_value( this->_flow_vars.v(), qp ) );
        if( this->_flow_vars.dim() == 3 )
          {
            U(2) = context.interior_value( this->_flow_vars.w(), qp );
          }

        // Compute the viscosity at this qp
        libMesh::Real _mu_qp = this->_mu(context, qp);

        libMesh::Real tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, this->_rho, U, _mu_qp, this->_is_steady );
        libMesh::Real tau_C = this->_stab_helper.compute_tau_continuity( tau_M, g );

        libMesh::RealGradient RM_s = this->_stab_helper.compute_res_momentum_steady( context, qp, this->_rho, _mu_qp );
        libMesh::Real RC = this->_stab_helper.compute_res_continuity( context, qp );

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) += ( - tau_C*RC*u_gradphi[i][qp](0)
                       - tau_M*RM_s(0)*this->_rho*U*u_gradphi[i][qp]  )*JxW[qp];

            Fv(i) += ( - tau_C*RC*u_gradphi[i][qp](1)
                       - tau_M*RM_s(1)*this->_rho*U*u_gradphi[i][qp] )*JxW[qp];

            if( this->_flow_vars.dim() == 3 )
              {
                (*Fw)(i) += ( - tau_C*RC*u_gradphi[i][qp](2)
                              - tau_M*RM_s(2)*this->_rho*U*u_gradphi[i][qp] )*JxW[qp];
              }
          }

        if( compute_jacobian )
          libmesh_not_implemented();
      }
  }

  template<class Mu>
  void IncompressibleNavierStokesSPGSMStabilization<Mu>::element_constraint
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.get_dof_indices(this->_press_var.p()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
      context.get_element_fe(this->_press_var.p())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(this->_press_var.p()); // R_{p}

    libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::RealGradient U( context.interior_value( this->_flow_vars.u(), qp ),
                                 context.interior_value( this->_flow_vars.v(), qp ) );
        if( this->_flow_vars.dim() == 3 )
          U(2) = context.interior_value( this->_flow_vars.w(), qp );

        // Compute the viscosity at this qp
        libMesh::Real _mu_qp = this->_mu(context, qp);

        libMesh::Real tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, this->_rho, U, _mu_qp, this->_is_steady );

        libMesh::RealGradient RM_s = this->_stab_helper.compute_res_momentum_steady( context, qp, this->_rho, _mu_qp );

        for (unsigned int i=0; i != n_p_dofs; i++)
          Fp(i) += tau_M*RM_s*p_dphi[i][qp]*JxW[qp];

        if( compute_jacobian )
          libmesh_not_implemented();
      }
  }

  template<class Mu>
  void IncompressibleNavierStokesSPGSMStabilization<Mu>::mass_residual
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.get_dof_indices(this->_press_var.p()).size();
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
      context.get_element_fe(this->_press_var.p())->get_dphi();

    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_flow_vars.u())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u()); // R_{p}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v()); // R_{p}
    libMesh::DenseSubVector<libMesh::Number> *Fw = NULL;


    if(this->_flow_vars.dim() == 3)
      Fw = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}

    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(this->_press_var.p()); // R_{p}

    libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::RealGradient U( context.fixed_interior_value( this->_flow_vars.u(), qp ),
                                 context.fixed_interior_value( this->_flow_vars.v(), qp ) );
        // Compute the viscosity at this qp
        libMesh::Real _mu_qp = this->_mu(context, qp);

        if( this->_flow_vars.dim() == 3 )
          {
            U(2) = context.fixed_interior_value( this->_flow_vars.w(), qp );
          }

        libMesh::Real tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, this->_rho, U, _mu_qp, false );

        libMesh::RealGradient RM_t = this->_stab_helper.compute_res_momentum_transient( context, qp, this->_rho );

        for (unsigned int i=0; i != n_p_dofs; i++)
          {
            Fp(i) -= tau_M*RM_t*p_dphi[i][qp]*JxW[qp];
          }

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) -= tau_M*RM_t(0)*this->_rho*U*u_gradphi[i][qp]*JxW[qp];

            Fv(i) -= tau_M*RM_t(1)*this->_rho*U*u_gradphi[i][qp]*JxW[qp];

            if( this->_flow_vars.dim() == 3 )
              {
                (*Fw)(i) -= tau_M*RM_t(2)*this->_rho*U*u_gradphi[i][qp]*JxW[qp];
              }
          }

        if( compute_jacobian )
          {
            libmesh_not_implemented();
          }

      }
  }

} // end namespace GRINS

// Instantiate
INSTANTIATE_INC_NS_SUBCLASS(IncompressibleNavierStokesSPGSMStabilization);
