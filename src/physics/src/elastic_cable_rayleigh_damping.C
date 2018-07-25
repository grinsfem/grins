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
#include "grins/elastic_cable_rayleigh_damping.h"

// GRINS
#include "grins/physics_naming.h"
#include "grins/elasticity_tensor.h"
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/quadrature.h"
#include "libmesh/first_order_unsteady_solver.h"

namespace GRINS
{
  template<typename StressStrainLaw>
  ElasticCableRayleighDamping<StressStrainLaw>::ElasticCableRayleighDamping( const PhysicsName& physics_name,
                                                                             const GetPot& input,
                                                                             bool is_compressible )
    : ElasticCableBase<StressStrainLaw>(physics_name,input,is_compressible),
    _lambda_factor(input("Physics/"+PhysicsNaming::elastic_cable_rayleigh_damping()+"/lambda_factor",0.0)),
    _mu_factor(input("Physics/"+PhysicsNaming::elastic_cable_rayleigh_damping()+"/mu_factor",0.0))
  {
    if( !input.have_variable("Physics/"+PhysicsNaming::elastic_cable_rayleigh_damping()+"/lambda_factor") )
      libmesh_error_msg("ERROR: Couldn't find Physics/"+PhysicsNaming::elastic_cable_rayleigh_damping()+"/lambda_factor in input!");

    if( !input.have_variable("Physics/"+PhysicsNaming::elastic_cable_rayleigh_damping()+"/mu_factor") )
      libmesh_error_msg("ERROR: Couldn't find Physics/"+PhysicsNaming::elastic_cable_rayleigh_damping()+"/mu_factor in input!");


    // If the user specified enabled subdomains in this Physics section,
    // that's an error; we're slave to ElasticCable.
    if( input.have_variable("Physics/"+PhysicsNaming::elastic_cable_rayleigh_damping()+"/enabled_subdomains" ) )
      libmesh_error_msg("ERROR: Cannot specify subdomains for "
                        +PhysicsNaming::elastic_cable_rayleigh_damping()
                        +"! Must specify subdomains through "
                        +PhysicsNaming::elastic_cable()+".");

    this->parse_enabled_subdomains(input,PhysicsNaming::elastic_cable());
  }

  template<typename StressStrainLaw>
  void ElasticCableRayleighDamping<StressStrainLaw>::auxiliary_init
  ( MultiphysicsSystem & system )
  {
    if( !this->is_steady() )
      {
        // Currently, we don't support first order time solvers
        const libMesh::TimeSolver & raw_time_solver = system.get_time_solver();

        const libMesh::FirstOrderUnsteadySolver * time_solver =
          dynamic_cast<const libMesh::FirstOrderUnsteadySolver *>(&raw_time_solver);

        if( time_solver )
          libmesh_error_msg("ERROR: First order time solvers not supported for ElasticCableRayleighDamping!");
      }
  }

  template<typename StressStrainLaw>
  void ElasticCableRayleighDamping<StressStrainLaw>::damping_residual
  ( bool compute_jacobian, AssemblyContext & context)
  {
    // First, do the "mass" contribution
    this->mass_residual_impl(compute_jacobian,
                             context,
                             &libMesh::FEMContext::interior_rate,
                             &libMesh::DiffContext::get_elem_solution_rate_derivative,
                             _mu_factor);

    // Now do the stiffness contribution
    const unsigned int n_u_dofs = context.get_dof_indices(this->_disp_vars.u()).size();

    const std::vector<libMesh::Real> &JxW =
      this->get_fe(context)->get_JxW();

    // Residuals, Jacobians that we're populating
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_disp_vars.u());
    libMesh::DenseSubVector<libMesh::Number>* Fv = NULL;
    libMesh::DenseSubVector<libMesh::Number>* Fw = NULL;

    libMesh::DenseSubMatrix<libMesh::Number>& Kuu = context.get_elem_jacobian(this->_disp_vars.u(),this->_disp_vars.u());
    libMesh::DenseSubMatrix<libMesh::Number>* Kuv = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kuw = NULL;

    libMesh::DenseSubMatrix<libMesh::Number>* Kvu = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kvv = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kvw = NULL;

    libMesh::DenseSubMatrix<libMesh::Number>* Kwu = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kwv = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kww = NULL;

    if( this->_disp_vars.dim() >= 2 )
      {
        Fv = &context.get_elem_residual(this->_disp_vars.v());

        Kuv = &context.get_elem_jacobian(this->_disp_vars.u(),this->_disp_vars.v());
        Kvu = &context.get_elem_jacobian(this->_disp_vars.v(),this->_disp_vars.u());
        Kvv = &context.get_elem_jacobian(this->_disp_vars.v(),this->_disp_vars.v());
      }

    if( this->_disp_vars.dim() == 3 )
      {
        Fw = &context.get_elem_residual(this->_disp_vars.w());

        Kuw = &context.get_elem_jacobian(this->_disp_vars.u(),this->_disp_vars.w());
        Kvw = &context.get_elem_jacobian(this->_disp_vars.v(),this->_disp_vars.w());
        Kwu = &context.get_elem_jacobian(this->_disp_vars.w(),this->_disp_vars.u());
        Kwv = &context.get_elem_jacobian(this->_disp_vars.w(),this->_disp_vars.v());
        Kww = &context.get_elem_jacobian(this->_disp_vars.w(),this->_disp_vars.w());
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    // All shape function gradients are w.r.t. master element coordinates
    const std::vector<std::vector<libMesh::Real> >& dphi_dxi = this->get_fe(context)->get_dphidxi();

    const libMesh::DenseSubVector<libMesh::Number>& u_coeffs = context.get_elem_solution( this->_disp_vars.u() );
    const libMesh::DenseSubVector<libMesh::Number>* v_coeffs = NULL;
    const libMesh::DenseSubVector<libMesh::Number>* w_coeffs = NULL;

    if( this->_disp_vars.dim() >= 2 )
      v_coeffs = &context.get_elem_solution( this->_disp_vars.v() );

    if( this->_disp_vars.dim() == 3 )
      w_coeffs = &context.get_elem_solution( this->_disp_vars.w() );

    const libMesh::DenseSubVector<libMesh::Number>& dudt_coeffs =
      context.get_elem_solution_rate( this->_disp_vars.u() );
    const libMesh::DenseSubVector<libMesh::Number>* dvdt_coeffs = NULL;
    const libMesh::DenseSubVector<libMesh::Number>* dwdt_coeffs = NULL;

    if( this->_disp_vars.dim() >= 2 )
      dvdt_coeffs = &context.get_elem_solution_rate( this->_disp_vars.v() );

    if( this->_disp_vars.dim() == 3 )
      dwdt_coeffs = &context.get_elem_solution_rate( this->_disp_vars.w() );

    // Need these to build up the covariant and contravariant metric tensors
    const std::vector<libMesh::RealGradient>& dxdxi  = this->get_fe(context)->get_dxyzdxi();

    const unsigned int dim = 1; // The cable dimension is always 1 for this physics

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // Gradients are w.r.t. master element coordinates
        libMesh::Gradient grad_u, grad_v, grad_w;
        libMesh::Gradient dgradu_dt, dgradv_dt, dgradw_dt;

        for( unsigned int d = 0; d < n_u_dofs; d++ )
          {
            libMesh::RealGradient u_gradphi( dphi_dxi[d][qp] );
            grad_u += u_coeffs(d)*u_gradphi;
            dgradu_dt += dudt_coeffs(d)*u_gradphi;

            if( this->_disp_vars.dim() >= 2 )
              {
                grad_v += (*v_coeffs)(d)*u_gradphi;
                dgradv_dt += (*dvdt_coeffs)(d)*u_gradphi;
              }

            if( this->_disp_vars.dim() == 3 )
              {
                grad_w += (*w_coeffs)(d)*u_gradphi;
                dgradw_dt += (*dwdt_coeffs)(d)*u_gradphi;
              }
          }

        libMesh::RealGradient grad_x( dxdxi[qp](0) );
        libMesh::RealGradient grad_y( dxdxi[qp](1) );
        libMesh::RealGradient grad_z( dxdxi[qp](2) );

        libMesh::TensorValue<libMesh::Real> a_cov, a_contra, A_cov, A_contra;
        libMesh::Real lambda_sq = 0;

        this->compute_metric_tensors( qp, *(this->get_fe(context)), context,
                                      grad_u, grad_v, grad_w,
                                      a_cov, a_contra, A_cov, A_contra,
                                      lambda_sq );

        // Compute stress tensor
        libMesh::TensorValue<libMesh::Real> tau;
        ElasticityTensor C;
        this->_stress_strain_law.compute_stress_and_elasticity(dim,a_contra,a_cov,A_contra,A_cov,tau,C);

        libMesh::Real jac = JxW[qp];

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            libMesh::RealGradient u_gradphi( dphi_dxi[i][qp] );

            libMesh::Real u_diag_factor = _lambda_factor*this->_A*jac*tau(0,0)*dgradu_dt(0)*u_gradphi(0);
            libMesh::Real v_diag_factor = _lambda_factor*this->_A*jac*tau(0,0)*dgradv_dt(0)*u_gradphi(0);
            libMesh::Real w_diag_factor = _lambda_factor*this->_A*jac*tau(0,0)*dgradw_dt(0)*u_gradphi(0);

            const libMesh::Real C1 = _lambda_factor*this->_A*jac*C(0,0,0,0)*u_gradphi(0);

            const libMesh::Real gamma_u = (grad_x(0)+grad_u(0));

            libMesh::Real gamma_v = 0.0;
            libMesh::Real gamma_w = 0.0;

            if( this->_disp_vars.dim() >= 2 )
              gamma_v = (grad_y(0)+grad_v(0));

            if( this->_disp_vars.dim() == 3 )
              gamma_w = (grad_z(0)+grad_w(0));

            const libMesh::Real x_term = C1*gamma_u;
            const libMesh::Real y_term = C1*gamma_v;
            const libMesh::Real z_term = C1*gamma_w;

            libMesh::Real dt_term = dgradu_dt(0)*gamma_u + dgradv_dt(0)*gamma_v + dgradw_dt(0)*gamma_w;

            Fu(i) += u_diag_factor + x_term*dt_term;

            if( this->_disp_vars.dim() >= 2 )
              (*Fv)(i) += v_diag_factor + y_term*dt_term;

            if( this->_disp_vars.dim() == 3 )
              (*Fw)(i) += w_diag_factor + z_term*dt_term;
          }

        if( compute_jacobian )
          {
            for(unsigned int i=0; i != n_u_dofs; i++)
              {
                libMesh::RealGradient u_gradphi_I( dphi_dxi[i][qp] );

                for(unsigned int j=0; j != n_u_dofs; j++)
                  {
                    libMesh::RealGradient u_gradphi_J( dphi_dxi[j][qp] );

                    libMesh::Real common_factor = _lambda_factor*this->_A*jac*u_gradphi_I(0);

                    const libMesh::Real diag_term_1 = common_factor*tau(0,0)*u_gradphi_J(0)*context.get_elem_solution_rate_derivative();

                    const libMesh::Real dgamma_du = ( u_gradphi_J(0)*(grad_x(0)+grad_u(0)) );

                    const libMesh::Real dgamma_dv = ( u_gradphi_J(0)*(grad_y(0)+grad_v(0)) );

                    const libMesh::Real dgamma_dw = ( u_gradphi_J(0)*(grad_z(0)+grad_w(0)) );

                    const libMesh::Real diag_term_2_factor = common_factor*C(0,0,0,0)*context.get_elem_solution_derivative();

                    Kuu(i,j) += diag_term_1 + dgradu_dt(0)*diag_term_2_factor*dgamma_du;

                    if( this->_disp_vars.dim() >= 2 )
                      {
                        (*Kuv)(i,j) += dgradu_dt(0)*diag_term_2_factor*dgamma_dv;
                        (*Kvu)(i,j) += dgradv_dt(0)*diag_term_2_factor*dgamma_du;
                        (*Kvv)(i,j) += diag_term_1 + dgradv_dt(0)*diag_term_2_factor*dgamma_dv;
                      }

                    if( this->_disp_vars.dim() == 3 )
                      {
                        (*Kuw)(i,j) += dgradu_dt(0)*diag_term_2_factor*dgamma_dw;
                        (*Kvw)(i,j) += dgradv_dt(0)*diag_term_2_factor*dgamma_dw;
                        (*Kwu)(i,j) += dgradw_dt(0)*diag_term_2_factor*dgamma_du;
                        (*Kwv)(i,j) += dgradw_dt(0)*diag_term_2_factor*dgamma_dv;
                        (*Kww)(i,j) += diag_term_1 + dgradw_dt(0)*diag_term_2_factor*dgamma_dw;
                      }

                    const libMesh::Real C1 = common_factor*C(0,0,0,0);

                    const libMesh::Real gamma_u = (grad_x(0)+grad_u(0));
                    const libMesh::Real gamma_v = (grad_y(0)+grad_v(0));
                    const libMesh::Real gamma_w = (grad_z(0)+grad_w(0));

                    const libMesh::Real x_term = C1*gamma_u;
                    const libMesh::Real y_term = C1*gamma_v;
                    const libMesh::Real z_term = C1*gamma_w;

                    const libMesh::Real ddtterm_du = u_gradphi_J(0)*(gamma_u*context.get_elem_solution_rate_derivative()
                                                                     + dgradu_dt(0)*context.get_elem_solution_derivative());

                    libMesh::Real ddtterm_dv = 0.0;
                    if( this->_disp_vars.dim() >= 2 )
                      ddtterm_dv = u_gradphi_J(0)*(gamma_v*context.get_elem_solution_rate_derivative()
                                                   + dgradv_dt(0)*context.get_elem_solution_derivative());

                    libMesh::Real ddtterm_dw = 0.0;
                    if( this->_disp_vars.dim() == 3 )
                      ddtterm_dw = u_gradphi_J(0)*(gamma_w*context.get_elem_solution_rate_derivative()
                                                   + dgradw_dt(0)*context.get_elem_solution_derivative());

                    Kuu(i,j) += x_term*ddtterm_du;

                    if( this->_disp_vars.dim() >= 2 )
                      {
                        (*Kuv)(i,j) += x_term*ddtterm_dv;
                        (*Kvu)(i,j) += y_term*ddtterm_du;
                        (*Kvv)(i,j) += y_term*ddtterm_dv;
                      }

                    if( this->_disp_vars.dim() == 3 )
                      {
                        (*Kuw)(i,j) += x_term*ddtterm_dw;
                        (*Kvw)(i,j) += y_term*ddtterm_dw;
                        (*Kwu)(i,j) += z_term*ddtterm_du;
                        (*Kwv)(i,j) += z_term*ddtterm_dv;
                        (*Kww)(i,j) += z_term*ddtterm_dw;
                      }

                    libMesh::Real dt_term = dgradu_dt(0)*gamma_u + dgradv_dt(0)*gamma_v + dgradw_dt(0)*gamma_w;

                    // Here, we're missing derivatives of C(0,0,0,0) w.r.t. strain
                    // Nonzero for hyperelasticity models
                    const libMesh::Real dxterm_du = C1*u_gradphi_J(0)*context.get_elem_solution_derivative();
                    const libMesh::Real dyterm_dv = dxterm_du;
                    const libMesh::Real dzterm_dw = dxterm_du;

                    Kuu(i,j) += dxterm_du*dt_term;

                    if( this->_disp_vars.dim() >= 2 )
                      (*Kvv)(i,j) += dyterm_dv*dt_term;

                    if( this->_disp_vars.dim() == 3 )
                      (*Kww)(i,j) += dzterm_dw*dt_term;

                  } // end j-loop
              } // end i-loop
          } // end if(compute_jacobian)
      } // end qp loop
  }

} // end namespace GRINS
