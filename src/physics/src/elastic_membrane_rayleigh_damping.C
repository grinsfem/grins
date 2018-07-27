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
#include "grins/elastic_membrane_rayleigh_damping.h"

// GRINS
#include "grins/physics_naming.h"
#include "grins/elasticity_tensor.h"
#include "grins/assembly_context.h"
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/quadrature.h"
#include "libmesh/first_order_unsteady_solver.h"

namespace GRINS
{
  template<typename StressStrainLaw>
  ElasticMembraneRayleighDamping<StressStrainLaw>::ElasticMembraneRayleighDamping( const PhysicsName& physics_name,
                                                                                   const GetPot& input,
                                                                                   bool is_compressible )
    : ElasticMembraneBase<StressStrainLaw>(physics_name,input,is_compressible),
    _lambda_factor(input("Physics/"+PhysicsNaming::elastic_membrane_rayleigh_damping()+"/lambda_factor",0.0)),
    _mu_factor(input("Physics/"+PhysicsNaming::elastic_membrane_rayleigh_damping()+"/mu_factor",0.0))
  {
    if( !input.have_variable("Physics/"+PhysicsNaming::elastic_membrane_rayleigh_damping()+"/lambda_factor") )
      libmesh_error_msg("ERROR: Couldn't find Physics/"+PhysicsNaming::elastic_membrane_rayleigh_damping()+"/lambda_factor in input!");

    if( !input.have_variable("Physics/"+PhysicsNaming::elastic_membrane_rayleigh_damping()+"/mu_factor") )
      libmesh_error_msg("ERROR: Couldn't find Physics/"+PhysicsNaming::elastic_membrane_rayleigh_damping()+"/mu_factor in input!");

    // If the user specified enabled subdomains in this Physics section,
    // that's an error; we're slave to ElasticMembrane.
    if( input.have_variable("Physics/"+PhysicsNaming::elastic_membrane_rayleigh_damping()+"/enabled_subdomains" ) )
      libmesh_error_msg("ERROR: Cannot specify subdomains for "
                        +PhysicsNaming::elastic_membrane_rayleigh_damping()
                        +"! Must specify subdomains through "
                        +PhysicsNaming::elastic_membrane()+".");

    this->parse_enabled_subdomains(input,PhysicsNaming::elastic_membrane());
  }

  template<typename StressStrainLaw>
  void ElasticMembraneRayleighDamping<StressStrainLaw>::auxiliary_init
  ( MultiphysicsSystem & system )
  {
    if( !this->is_steady() )
      {
        // Currently, we don't support first order time solvers
        const libMesh::TimeSolver & raw_time_solver = system.get_time_solver();

        const libMesh::FirstOrderUnsteadySolver * time_solver =
          dynamic_cast<const libMesh::FirstOrderUnsteadySolver *>(&raw_time_solver);

        if( time_solver )
          libmesh_error_msg("ERROR: First order time solvers not supported for ElasticMembraneRayleighDamping!");
      }
  }

  template<typename StressStrainLaw>
  void ElasticMembraneRayleighDamping<StressStrainLaw>::damping_residual
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // First do the mass part
    this->mass_residual_impl(compute_jacobian,
                             context,
                             &libMesh::FEMContext::interior_rate,
                             &libMesh::DiffContext::get_elem_solution_rate_derivative,
                             _mu_factor);

    // Now do the stiffness part
    const unsigned int n_u_dofs = context.get_dof_indices(this->_disp_vars.u()).size();

    const std::vector<libMesh::Real> &JxW =
      this->get_fe(context)->get_JxW();

    // Residuals that we're populating
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_disp_vars.u());
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_disp_vars.v());
    libMesh::DenseSubVector<libMesh::Number>* Fw = NULL;

    libMesh::DenseSubMatrix<libMesh::Number>& Kuu = context.get_elem_jacobian(this->_disp_vars.u(),this->_disp_vars.u());
    libMesh::DenseSubMatrix<libMesh::Number>& Kuv = context.get_elem_jacobian(this->_disp_vars.u(),this->_disp_vars.v());

    libMesh::DenseSubMatrix<libMesh::Number>& Kvu = context.get_elem_jacobian(this->_disp_vars.v(),this->_disp_vars.u());
    libMesh::DenseSubMatrix<libMesh::Number>& Kvv = context.get_elem_jacobian(this->_disp_vars.v(),this->_disp_vars.v());

    libMesh::DenseSubMatrix<libMesh::Number>* Kuw = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kvw = NULL;

    libMesh::DenseSubMatrix<libMesh::Number>* Kwu = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kwv = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kww = NULL;

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
    const std::vector<std::vector<libMesh::Real> >& dphi_dxi =
      this->get_fe(context)->get_dphidxi();

    const std::vector<std::vector<libMesh::Real> >& dphi_deta =
      this->get_fe(context)->get_dphideta();

    const libMesh::DenseSubVector<libMesh::Number>& u_coeffs = context.get_elem_solution( this->_disp_vars.u() );
    const libMesh::DenseSubVector<libMesh::Number>& v_coeffs = context.get_elem_solution( this->_disp_vars.v() );
    const libMesh::DenseSubVector<libMesh::Number>* w_coeffs = NULL;

    const libMesh::DenseSubVector<libMesh::Number>& dudt_coeffs = context.get_elem_solution_rate( this->_disp_vars.u() );
    const libMesh::DenseSubVector<libMesh::Number>& dvdt_coeffs = context.get_elem_solution_rate( this->_disp_vars.v() );
    const libMesh::DenseSubVector<libMesh::Number>* dwdt_coeffs = NULL;

    if( this->_disp_vars.dim() == 3 )
      {
        w_coeffs = &context.get_elem_solution( this->_disp_vars.w() );
        dwdt_coeffs = &context.get_elem_solution_rate( this->_disp_vars.w() );
      }

    // Need these to build up the covariant and contravariant metric tensors
    const std::vector<libMesh::RealGradient>& dxdxi  = this->get_fe(context)->get_dxyzdxi();
    const std::vector<libMesh::RealGradient>& dxdeta = this->get_fe(context)->get_dxyzdeta();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // Gradients are w.r.t. master element coordinates
        libMesh::Gradient grad_u, grad_v, grad_w;
        libMesh::Gradient dgradu_dt, dgradv_dt, dgradw_dt;

        for( unsigned int d = 0; d < n_u_dofs; d++ )
          {
            libMesh::RealGradient u_gradphi( dphi_dxi[d][qp], dphi_deta[d][qp] );
            grad_u += u_coeffs(d)*u_gradphi;
            grad_v += v_coeffs(d)*u_gradphi;

            dgradu_dt += dudt_coeffs(d)*u_gradphi;
            dgradv_dt += dvdt_coeffs(d)*u_gradphi;

            if( this->_disp_vars.dim() == 3 )
              {
                grad_w += (*w_coeffs)(d)*u_gradphi;
                dgradw_dt += (*dwdt_coeffs)(d)*u_gradphi;
              }
          }

        libMesh::RealGradient grad_x( dxdxi[qp](0), dxdeta[qp](0) );
        libMesh::RealGradient grad_y( dxdxi[qp](1), dxdeta[qp](1) );
        libMesh::RealGradient grad_z( dxdxi[qp](2), dxdeta[qp](2) );


        libMesh::TensorValue<libMesh::Real> a_cov, a_contra, A_cov, A_contra;
        libMesh::Real lambda_sq = 0;

        this->compute_metric_tensors( qp, *(this->get_fe(context)), context,
                                      grad_u, grad_v, grad_w,
                                      a_cov, a_contra, A_cov, A_contra,
                                      lambda_sq );

        const unsigned int manifold_dim = 2; // The manifold dimension is always 2 for this physics

        // Compute stress and elasticity tensors
        libMesh::TensorValue<libMesh::Real> tau;
        ElasticityTensor C;
        this->_stress_strain_law.compute_stress_and_elasticity(manifold_dim,a_contra,a_cov,A_contra,A_cov,tau,C);

        libMesh::Real jac = JxW[qp];

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            libMesh::RealGradient u_gradphi( dphi_dxi[i][qp], dphi_deta[i][qp] );

            for( unsigned int alpha = 0; alpha < manifold_dim; alpha++ )
              {
                for( unsigned int beta = 0; beta < manifold_dim; beta++ )
                  {
                    const libMesh::Real common_factor = _lambda_factor*0.5*this->_h0*jac;
                    const libMesh::Real factor = common_factor*tau(alpha,beta);

                    libMesh::Real u_diagterm = dgradu_dt(beta)*u_gradphi(alpha) + dgradu_dt(alpha)*u_gradphi(beta);
                    libMesh::Real v_diagterm = dgradv_dt(beta)*u_gradphi(alpha) + dgradv_dt(alpha)*u_gradphi(beta);
                    libMesh::Real w_diagterm = dgradw_dt(beta)*u_gradphi(alpha) + dgradw_dt(alpha)*u_gradphi(beta);

                    Fu(i) += factor*u_diagterm;
                    Fv(i) += factor*v_diagterm;

                    if( this->_disp_vars.dim() == 3 )
                      (*Fw)(i) += factor*w_diagterm;

                    for( unsigned int lambda = 0; lambda < manifold_dim; lambda++ )
                      {
                        for( unsigned int mu = 0; mu < manifold_dim; mu++ )
                          {
                            const libMesh::Real C1 = common_factor*C(alpha,beta,lambda,mu);

                            const libMesh::Real gamma_u = 0.5*( dgradu_dt(lambda)*(grad_x(mu)+grad_u(mu)) +
                                                                (grad_x(lambda)+grad_u(lambda))*dgradu_dt(mu) );

                            const libMesh::Real gamma_v = 0.5*( dgradv_dt(lambda)*(grad_y(mu)+grad_v(mu)) +
                                                                (grad_y(lambda)+grad_v(lambda))*dgradv_dt(mu) );

                            const libMesh::Real x_term = C1*( (grad_x(beta)+grad_u(beta))*u_gradphi(alpha) +
                                                              (grad_x(alpha)+grad_u(alpha))*u_gradphi(beta) );

                            const libMesh::Real y_term = C1*( (grad_y(beta)+grad_v(beta))*u_gradphi(alpha) +
                                                              (grad_y(alpha)+grad_v(alpha))*u_gradphi(beta) );

                            libMesh::Real gamma_sum = gamma_u + gamma_v;
                            if( this->_disp_vars.dim() == 3 )
                              {
                                const libMesh::Real gamma_w = 0.5*( dgradw_dt(lambda)*(grad_z(mu)+grad_w(mu)) +
                                                                    (grad_z(lambda)+grad_w(lambda))*dgradw_dt(mu) );
                                gamma_sum += gamma_w;
                              }

                            Fu(i) += x_term*gamma_sum;

                            Fv(i) += y_term*gamma_sum;

                            if( this->_disp_vars.dim() == 3 )
                              {
                                const libMesh::Real z_term = C1*( (grad_z(beta)+grad_w(beta))*u_gradphi(alpha) +
                                                                  (grad_z(alpha)+grad_w(alpha))*u_gradphi(beta) );

                                (*Fw)(i) += z_term*gamma_sum;
                              }
                          }
                      }
                  }
              }
          }

        if( compute_jacobian )
          {
            for (unsigned int i=0; i != n_u_dofs; i++)
              {
                libMesh::RealGradient u_gradphi_i( dphi_dxi[i][qp], dphi_deta[i][qp] );

                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    libMesh::RealGradient u_gradphi_j( dphi_dxi[j][qp], dphi_deta[j][qp] );

                    for( unsigned int alpha = 0; alpha < manifold_dim; alpha++ )
                      {
                        for( unsigned int beta = 0; beta < manifold_dim; beta++ )
                          {
                            const libMesh::Real common_factor = _lambda_factor*0.5*this->_h0*jac;

                            const libMesh::Real factor = common_factor*tau(alpha,beta);

                            const libMesh::Real ddiagterm = u_gradphi_j(beta)*u_gradphi_i(alpha) + u_gradphi_j(alpha)*u_gradphi_i(beta);

                            Kuu(i,j) += factor*ddiagterm*context.get_elem_solution_rate_derivative();
                            Kvv(i,j) += factor*ddiagterm*context.get_elem_solution_rate_derivative();

                            if( this->_disp_vars.dim() == 3 )
                              (*Kww)(i,j) += factor*ddiagterm*context.get_elem_solution_rate_derivative();

                            libMesh::Real u_diagterm = dgradu_dt(beta)*u_gradphi_i(alpha) + dgradu_dt(alpha)*u_gradphi_i(beta);
                            libMesh::Real v_diagterm = dgradv_dt(beta)*u_gradphi_i(alpha) + dgradv_dt(alpha)*u_gradphi_i(beta);
                            libMesh::Real w_diagterm = 0.0;
                            if( this->_disp_vars.dim() == 3 )
                              w_diagterm = dgradw_dt(beta)*u_gradphi_i(alpha) + dgradw_dt(alpha)*u_gradphi_i(beta);

                            for( unsigned int lambda = 0; lambda < manifold_dim; lambda++ )
                              {
                                for( unsigned int mu = 0; mu < manifold_dim; mu++ )
                                  {
                                    const libMesh::Real dgamma_du = 0.5*( u_gradphi_j(lambda)*(grad_x(mu)+grad_u(mu)) +
                                                                          (grad_x(lambda)+grad_u(lambda))*u_gradphi_j(mu) );

                                    const libMesh::Real dgamma_dv = 0.5*( u_gradphi_j(lambda)*(grad_y(mu)+grad_v(mu)) +
                                                                          (grad_y(lambda)+grad_v(lambda))*u_gradphi_j(mu) );

                                    const libMesh::Real C1 = common_factor*C(alpha,beta,lambda,mu);

                                    const libMesh::Real dfactor_du = C1*dgamma_du*context.get_elem_solution_derivative();
                                    const libMesh::Real dfactor_dv = C1*dgamma_dv*context.get_elem_solution_derivative();

                                    Kuu(i,j) += u_diagterm*dfactor_du;
                                    Kuv(i,j) += u_diagterm*dfactor_dv;

                                    Kvu(i,j) += v_diagterm*dfactor_du;
                                    Kvv(i,j) += v_diagterm*dfactor_dv;

                                    if( this->_disp_vars.dim() == 3 )
                                      {
                                        const libMesh::Real dgamma_dw = 0.5*( u_gradphi_j(lambda)*(grad_z(mu)+grad_w(mu)) +
                                                                              (grad_z(lambda)+grad_w(lambda))*u_gradphi_j(mu) );

                                        const libMesh::Real dfactor_dw = C1*dgamma_dw*context.get_elem_solution_derivative();

                                        (*Kuw)(i,j) += u_diagterm*dfactor_dw;

                                        (*Kvw)(i,j) += v_diagterm*dfactor_dw;

                                        (*Kwu)(i,j) += w_diagterm*dfactor_du;
                                        (*Kwv)(i,j) += w_diagterm*dfactor_dv;
                                        (*Kww)(i,j) += w_diagterm*dfactor_dw;
                                      }

                                    const libMesh::Real gamma_u = 0.5*( dgradu_dt(lambda)*(grad_x(mu)+grad_u(mu)) +
                                                                        (grad_x(lambda)+grad_u(lambda))*dgradu_dt(mu) );

                                    const libMesh::Real gamma_v = 0.5*( dgradv_dt(lambda)*(grad_y(mu)+grad_v(mu)) +
                                                                        (grad_y(lambda)+grad_v(lambda))*dgradv_dt(mu) );

                                    const libMesh::Real x_term = C1*( (grad_x(beta)+grad_u(beta))*u_gradphi_i(alpha) +
                                                                      (grad_x(alpha)+grad_u(alpha))*u_gradphi_i(beta) );

                                    const libMesh::Real y_term = C1*( (grad_y(beta)+grad_v(beta))*u_gradphi_i(alpha) +
                                                                      (grad_y(alpha)+grad_v(alpha))*u_gradphi_i(beta) );

                                    libMesh::Real gamma_sum = gamma_u + gamma_v;
                                    if( this->_disp_vars.dim() == 3 )
                                      {
                                        const libMesh::Real gamma_w = 0.5*( dgradw_dt(lambda)*(grad_z(mu)+grad_w(mu)) +
                                                                            (grad_z(lambda)+grad_w(lambda))*dgradw_dt(mu) );
                                        gamma_sum += gamma_w;
                                      }

                                    // Here, we're missing derivatives of C(alpha,beta,lambda,mu) w.r.t. strain
                                    // Nonzero for hyperelasticity models
                                    const libMesh::Real dxterm_du = C1*(  u_gradphi_j(beta)*u_gradphi_i(alpha)
                                                                          + u_gradphi_j(alpha)*u_gradphi_i(beta) )*context.get_elem_solution_derivative();

                                    const libMesh::Real dyterm_dv = dxterm_du;

                                    Kuu(i,j) += gamma_sum*dxterm_du;

                                    Kvv(i,j) += gamma_sum*dyterm_dv;

                                    libMesh::Real z_term = 0.0;

                                    if( this->_disp_vars.dim() == 3 )
                                      {
                                        z_term = C1*( (grad_z(beta)+grad_w(beta))*u_gradphi_i(alpha) +
                                                      (grad_z(alpha)+grad_w(alpha))*u_gradphi_i(beta) );

                                        // Here, we're missing derivatives of C(alpha,beta,lambda,mu) w.r.t. strain
                                        // Nonzero for hyperelasticity models
                                        const libMesh::Real dzterm_dw = dxterm_du;

                                        (*Kww)(i,j) += gamma_sum*dzterm_dw;
                                      }

                                    const libMesh::Real dgamma_sum_du = 0.5*(   u_gradphi_j(lambda)*(grad_x(mu)+grad_u(mu))*context.get_elem_solution_rate_derivative()
                                                                                + dgradu_dt(lambda)*u_gradphi_j(mu)*context.get_elem_solution_derivative()
                                                                                + u_gradphi_j(mu)*(grad_x(lambda)+grad_u(lambda))*context.get_elem_solution_rate_derivative()
                                                                                + dgradu_dt(mu)*u_gradphi_j(lambda)*context.get_elem_solution_derivative() );

                                    const libMesh::Real dgamma_sum_dv = 0.5*(   u_gradphi_j(lambda)*(grad_y(mu)+grad_v(mu))*context.get_elem_solution_rate_derivative()
                                                                                + dgradv_dt(lambda)*u_gradphi_j(mu)*context.get_elem_solution_derivative()
                                                                                + u_gradphi_j(mu)*(grad_y(lambda)+grad_v(lambda))*context.get_elem_solution_rate_derivative()
                                                                                + dgradv_dt(mu)*u_gradphi_j(lambda)*context.get_elem_solution_derivative() );

                                    Kuu(i,j) += x_term*dgamma_sum_du;
                                    Kuv(i,j) += x_term*dgamma_sum_dv;

                                    Kvu(i,j) += y_term*dgamma_sum_du;
                                    Kvv(i,j) += y_term*dgamma_sum_dv;

                                    if( this->_disp_vars.dim() == 3 )
                                      {
                                        const libMesh::Real dgamma_sum_dw = 0.5*(   u_gradphi_j(lambda)*(grad_z(mu)+grad_w(mu))*context.get_elem_solution_rate_derivative()
                                                                                    + dgradw_dt(lambda)*u_gradphi_j(mu)*context.get_elem_solution_derivative()
                                                                                    + u_gradphi_j(mu)*(grad_z(lambda)+grad_w(lambda))*context.get_elem_solution_rate_derivative()
                                                                                    + dgradw_dt(mu)*u_gradphi_j(lambda)*context.get_elem_solution_derivative() );

                                        (*Kuw)(i,j) += x_term*dgamma_sum_dw;
                                        (*Kvw)(i,j) += y_term*dgamma_sum_dw;
                                        (*Kwu)(i,j) += z_term*dgamma_sum_du;
                                        (*Kwv)(i,j) += z_term*dgamma_sum_dv;
                                        (*Kww)(i,j) += z_term*dgamma_sum_dw;
                                      }
                                  }
                              }
                          }
                      }
                  }
              }
          }

      }
  }


} // end namespace GRINS
