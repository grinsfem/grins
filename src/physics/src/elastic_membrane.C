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
#include "grins/elastic_membrane.h"

// GRINS
#include "grins_config.h"
#include "grins/math_constants.h"
#include "grins/assembly_context.h"
#include "grins/generic_ic_handler.h"
#include "grins/postprocessed_quantities.h"
#include "grins/elasticity_tensor.h"
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/quadrature.h"
#include "libmesh/boundary_info.h"
#include "libmesh/fem_system.h"
#include "libmesh/elem.h"

namespace GRINS
{
  template<typename StressStrainLaw>
  ElasticMembrane<StressStrainLaw>::ElasticMembrane( const GRINS::PhysicsName& physics_name, const GetPot& input,
                                                     bool is_compressible )
    : ElasticMembraneBase<StressStrainLaw>(physics_name,input,is_compressible)
  {
    this->_ic_handler = new GenericICHandler(physics_name, input);
  }

  template<typename StressStrainLaw>
  void ElasticMembrane<StressStrainLaw>::register_postprocessing_vars( const GetPot& input,
                                                                       PostProcessedQuantities<libMesh::Real>& postprocessing )
  {
    std::string section = "Physics/"+PhysicsNaming::elastic_membrane()+"/output_vars";

    if( input.have_variable(section) )
      {
        unsigned int n_vars = input.vector_variable_size(section);

        for( unsigned int v = 0; v < n_vars; v++ )
          {
            std::string name = input(section,"DIE!",v);

            if( name == std::string("stress") )
              {
                // sigma_xx, sigma_xy, sigma_yy, sigma_yx = sigma_xy
                // sigma_zz = 0 by assumption of this Physics
                _stress_indices.resize(7);

                this->_stress_indices[0] = postprocessing.register_quantity("stress_xx");

                this->_stress_indices[1] = postprocessing.register_quantity("stress_xy");

                this->_stress_indices[2] = postprocessing.register_quantity("stress_yy");

                this->_stress_indices[3] = postprocessing.register_quantity("sigma_max");

                this->_stress_indices[4] = postprocessing.register_quantity("sigma_1");

                this->_stress_indices[5] = postprocessing.register_quantity("sigma_2");

                this->_stress_indices[6] = postprocessing.register_quantity("sigma_3");
              }
            else if( name == std::string( "stress_zz" ) )
              {
                // This is mostly for sanity checking the plane stress condition
                this->_stress_zz_index = postprocessing.register_quantity("stress_zz");
              }
            else if( name == std::string("strain") )
              {
                // eps_xx, eps_xy, eps_yy, eps_yx = eps_xy
                _strain_indices.resize(3);

                this->_strain_indices[0] = postprocessing.register_quantity("strain_xx");

                this->_strain_indices[1] = postprocessing.register_quantity("strain_xy");

                this->_strain_indices[2] = postprocessing.register_quantity("strain_yy");
              }
            else
              {
                std::cerr << "Error: Invalue output_vars value for "+PhysicsNaming::elastic_membrane() << std::endl
                          << "       Found " << name << std::endl
                          << "       Acceptable values are: stress" << std::endl
                          << "                              strain" << std::endl;
                libmesh_error();
              }
          }
      }
  }

  template<typename StressStrainLaw>
  void ElasticMembrane<StressStrainLaw>::element_time_derivative
  ( bool compute_jacobian, AssemblyContext & context )
  {
    unsigned int u_var = this->_disp_vars.u();
    unsigned int v_var = this->_disp_vars.v();

    const unsigned int n_u_dofs = context.get_dof_indices(u_var).size();

    const std::vector<libMesh::Real> &JxW =
      this->get_fe(context)->get_JxW();

    const MultiphysicsSystem & system = context.get_multiphysics_system();

    unsigned int u_dot_var = system.get_second_order_dot_var(u_var);
    unsigned int v_dot_var = system.get_second_order_dot_var(v_var);

    // Residuals that we're populating
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(u_dot_var);
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(v_dot_var);
    libMesh::DenseSubVector<libMesh::Number>* Fw = NULL;

    libMesh::DenseSubMatrix<libMesh::Number>& Kuu = context.get_elem_jacobian(u_dot_var,u_var);
    libMesh::DenseSubMatrix<libMesh::Number>& Kuv = context.get_elem_jacobian(u_dot_var,v_var);

    libMesh::DenseSubMatrix<libMesh::Number>& Kvu = context.get_elem_jacobian(v_dot_var,u_var);
    libMesh::DenseSubMatrix<libMesh::Number>& Kvv = context.get_elem_jacobian(v_dot_var,v_var);

    libMesh::DenseSubMatrix<libMesh::Number>* Kuw = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kvw = NULL;

    libMesh::DenseSubMatrix<libMesh::Number>* Kwu = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kwv = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kww = NULL;

    unsigned int w_var = libMesh::invalid_uint;
    unsigned int w_dot_var = libMesh::invalid_uint;
    if( this->_disp_vars.dim() == 3 )
      {
        w_var = this->_disp_vars.w();
        w_dot_var = system.get_second_order_dot_var(w_var);

        Fw = &context.get_elem_residual(w_dot_var);

        Kuw = &context.get_elem_jacobian(u_dot_var,w_var);
        Kvw = &context.get_elem_jacobian(v_dot_var,w_var);
        Kwu = &context.get_elem_jacobian(w_dot_var,u_var);
        Kwv = &context.get_elem_jacobian(w_dot_var,v_var);
        Kww = &context.get_elem_jacobian(w_dot_var,w_var);
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    // All shape function gradients are w.r.t. master element coordinates
    const std::vector<std::vector<libMesh::Real> >& dphi_dxi =
      this->get_fe(context)->get_dphidxi();

    const std::vector<std::vector<libMesh::Real> >& dphi_deta =
      this->get_fe(context)->get_dphideta();

    // Need these to build up the covariant and contravariant metric tensors
    const std::vector<libMesh::RealGradient>& dxdxi  = this->get_fe(context)->get_dxyzdxi();
    const std::vector<libMesh::RealGradient>& dxdeta = this->get_fe(context)->get_dxyzdeta();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Gradient grad_u, grad_v,grad_w;
        this->get_grad_disp(context, qp, grad_u,grad_v,grad_w);

        // Compute stress and elasticity tensors
        libMesh::TensorValue<libMesh::Real> tau;
        ElasticityTensor C;
        this->get_stress_and_elasticity(context,qp,grad_u,grad_v,grad_w,tau,C);

        libMesh::RealGradient grad_x( dxdxi[qp](0), dxdeta[qp](0) );
        libMesh::RealGradient grad_y( dxdxi[qp](1), dxdeta[qp](1) );
        libMesh::RealGradient grad_z( dxdxi[qp](2), dxdeta[qp](2) );

        const unsigned int manifold_dim = 2; // The manifold dimension is always 2 for this physics

        libMesh::Real jac = JxW[qp];

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            libMesh::RealGradient u_gradphi( dphi_dxi[i][qp], dphi_deta[i][qp] );

            for( unsigned int alpha = 0; alpha < manifold_dim; alpha++ )
              {
                for( unsigned int beta = 0; beta < manifold_dim; beta++ )
                  {
                    libMesh::Real factor = 0.5*tau(alpha,beta)*this->_h0*jac;

                    Fu(i) += factor*( (grad_x(beta) + grad_u(beta))*u_gradphi(alpha) +
                                      (grad_x(alpha) + grad_u(alpha))*u_gradphi(beta) );

                    Fv(i) += factor*( (grad_y(beta) + grad_v(beta))*u_gradphi(alpha) +
                                      (grad_y(alpha) + grad_v(alpha))*u_gradphi(beta) );

                    if( this->_disp_vars.dim() == 3 )
                      (*Fw)(i) += factor*( (grad_z(beta) + grad_w(beta))*u_gradphi(alpha) +
                                           (grad_z(alpha) + grad_w(alpha))*u_gradphi(beta) );
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
                            const libMesh::Real diag_term = 0.5*this->_h0*jac*tau(alpha,beta)*context.get_elem_solution_derivative()*
                              ( u_gradphi_j(beta)*u_gradphi_i(alpha) +
                                u_gradphi_j(alpha)*u_gradphi_i(beta) );
                            Kuu(i,j) += diag_term;

                            Kvv(i,j) += diag_term;

                            if( this->_disp_vars.dim() == 3 )
                              (*Kww)(i,j) += diag_term;

                            for( unsigned int lambda = 0; lambda < manifold_dim; lambda++ )
                              {
                                for( unsigned int mu = 0; mu < manifold_dim; mu++ )
                                  {
                                    const libMesh::Real dgamma_du = 0.5*( u_gradphi_j(lambda)*(grad_x(mu)+grad_u(mu)) +
                                                                          (grad_x(lambda)+grad_u(lambda))*u_gradphi_j(mu) );

                                    const libMesh::Real dgamma_dv = 0.5*( u_gradphi_j(lambda)*(grad_y(mu)+grad_v(mu)) +
                                                                          (grad_y(lambda)+grad_v(lambda))*u_gradphi_j(mu) );

                                    const libMesh::Real C1 = 0.5*this->_h0*jac*C(alpha,beta,lambda,mu)*context.get_elem_solution_derivative();

                                    const libMesh::Real x_term = C1*( (grad_x(beta)+grad_u(beta))*u_gradphi_i(alpha) +
                                                                      (grad_x(alpha)+grad_u(alpha))*u_gradphi_i(beta) );

                                    const libMesh::Real y_term = C1*( (grad_y(beta)+grad_v(beta))*u_gradphi_i(alpha) +
                                                                      (grad_y(alpha)+grad_v(alpha))*u_gradphi_i(beta) );

                                    Kuu(i,j) += x_term*dgamma_du;

                                    Kuv(i,j) += x_term*dgamma_dv;

                                    Kvu(i,j) += y_term*dgamma_du;

                                    Kvv(i,j) += y_term*dgamma_dv;

                                    if( this->_disp_vars.dim() == 3 )
                                      {
                                        const libMesh::Real dgamma_dw = 0.5*( u_gradphi_j(lambda)*(grad_z(mu)+grad_w(mu)) +
                                                                              (grad_z(lambda)+grad_w(lambda))*u_gradphi_j(mu) );

                                        const libMesh::Real z_term = C1*( (grad_z(beta)+grad_w(beta))*u_gradphi_i(alpha) +
                                                                          (grad_z(alpha)+grad_w(alpha))*u_gradphi_i(beta) );

                                        (*Kuw)(i,j) += x_term*dgamma_dw;

                                        (*Kvw)(i,j) += y_term*dgamma_dw;

                                        (*Kwu)(i,j) += z_term*dgamma_du;

                                        (*Kwv)(i,j) += z_term*dgamma_dv;

                                        (*Kww)(i,j) += z_term*dgamma_dw;
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

  template<typename StressStrainLaw>
  void ElasticMembrane<StressStrainLaw>::element_constraint
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // Only compute the constraint is tracking lambda_sq as an independent variable
    if( this->_is_compressible )
      {
        unsigned int n_qpoints = context.get_element_qrule().n_points();

        const std::vector<libMesh::Real> & JxW = context.get_element_fe(this->_lambda_sq_var)->get_JxW();

        libMesh::DenseSubVector<libMesh::Number> & Fl = context.get_elem_residual(this->_lambda_sq_var);

        const std::vector<std::vector<libMesh::Real> >& phi =
          context.get_element_fe(this->_lambda_sq_var)->get_phi();

        const unsigned int n_lambda_sq_dofs = context.get_dof_indices(this->_lambda_sq_var).size();

        for (unsigned int qp=0; qp != n_qpoints; qp++)
          {
            libMesh::Real jac = JxW[qp];

            libMesh::Gradient grad_u, grad_v,grad_w;
            this->get_grad_disp(context, qp, grad_u,grad_v,grad_w);

            libMesh::TensorValue<libMesh::Real> a_cov, a_contra, A_cov, A_contra;
            libMesh::Real lambda_sq = 0;

            this->compute_metric_tensors( qp, *(this->get_fe(context)), context,
                                          grad_u, grad_v, grad_w,
                                          a_cov, a_contra, A_cov, A_contra,
                                          lambda_sq );

            libMesh::Real stress_33 = this->_stress_strain_law.compute_33_stress( a_contra, a_cov, A_contra, A_cov );

            for (unsigned int i=0; i != n_lambda_sq_dofs; i++)
              Fl(i) += stress_33*phi[i][qp]*jac;

            if( compute_jacobian )
              libmesh_not_implemented();
          }
      } // is_compressible
  }

  template<typename StressStrainLaw>
  void ElasticMembrane<StressStrainLaw>::compute_postprocessed_quantity( unsigned int quantity_index,
                                                                         const AssemblyContext& context,
                                                                         const libMesh::Point& point,
                                                                         libMesh::Real& value )
  {
    bool is_stress = false;
    if( !_stress_indices.empty() )
      is_stress= ( _stress_indices[0] == quantity_index ||
                   _stress_indices[1] == quantity_index ||
                   _stress_indices[2] == quantity_index ||
                   _stress_indices[3] == quantity_index ||
                   _stress_indices[4] == quantity_index ||
                   _stress_indices[5] == quantity_index ||
                   _stress_indices[6] == quantity_index ||
                   _stress_zz_index == quantity_index   );

    bool is_strain = false;
    if( !_strain_indices.empty() )
      is_strain = ( _strain_indices[0] == quantity_index ||
                    _strain_indices[1] == quantity_index ||
                    _strain_indices[2] == quantity_index   );

    if( is_stress || is_strain )
      {
        const unsigned int n_u_dofs = context.get_dof_indices(this->_disp_vars.u()).size();

        const libMesh::DenseSubVector<libMesh::Number>& u_coeffs = context.get_elem_solution( this->_disp_vars.u() );
        const libMesh::DenseSubVector<libMesh::Number>& v_coeffs = context.get_elem_solution( this->_disp_vars.v() );
        const libMesh::DenseSubVector<libMesh::Number>* w_coeffs = NULL;

        if( this->_disp_vars.dim() == 3 )
          w_coeffs = &context.get_elem_solution( this->_disp_vars.w() );

        // Build new FE for the current point. We need this to build tensors at point.
        std::unique_ptr<libMesh::FEGenericBase<libMesh::Real> > fe_new =
          this->build_new_fe( &context.get_elem(), this->get_fe(context),
                              point );

        const std::vector<std::vector<libMesh::Real> >& dphi_dxi =
          fe_new->get_dphidxi();

        const std::vector<std::vector<libMesh::Real> >& dphi_deta =
          fe_new->get_dphideta();

        // Need these to build up the covariant and contravariant metric tensors
        const std::vector<libMesh::RealGradient>& dxdxi  = fe_new->get_dxyzdxi();
        const std::vector<libMesh::RealGradient>& dxdeta = fe_new->get_dxyzdeta();

        // Gradients are w.r.t. master element coordinates
        libMesh::Gradient grad_u, grad_v, grad_w;
        for( unsigned int d = 0; d < n_u_dofs; d++ )
          {
            libMesh::RealGradient u_gradphi( dphi_dxi[d][0], dphi_deta[d][0] );
            grad_u += u_coeffs(d)*u_gradphi;
            grad_v += v_coeffs(d)*u_gradphi;

            if( this->_disp_vars.dim() == 3 )
              grad_w += (*w_coeffs)(d)*u_gradphi;
          }

        libMesh::RealGradient grad_x( dxdxi[0](0), dxdeta[0](0) );
        libMesh::RealGradient grad_y( dxdxi[0](1), dxdeta[0](1) );
        libMesh::RealGradient grad_z( dxdxi[0](2), dxdeta[0](2) );

        libMesh::TensorValue<libMesh::Real> a_cov, a_contra, A_cov, A_contra;
        libMesh::Real lambda_sq = 0;

        // We're only computing one point at a time, so qp = 0 always
        this->compute_metric_tensors(0, *fe_new, context, grad_u, grad_v, grad_w,
                                     a_cov, a_contra, A_cov, A_contra, lambda_sq );

        // We have everything we need for strain now, so check if we are computing strain
        if( is_strain )
          {
            if( _strain_indices[0] == quantity_index )
              {
                value = 0.5*(A_cov(0,0) - a_cov(0,0));
              }
            else if( _strain_indices[1] == quantity_index )
              {
                value = 0.5*(A_cov(0,1) - a_cov(0,1));
              }
            else if( _strain_indices[2] == quantity_index )
              {
                value = 0.5*(A_cov(1,1) - a_cov(1,1));
              }
            else
              {
                //Wat?!
                libmesh_error();
              }
            return;
          }

        if( is_stress )
          {
            libMesh::Real det_a = a_cov(0,0)*a_cov(1,1) - a_cov(0,1)*a_cov(1,0);
            libMesh::Real det_A = A_cov(0,0)*A_cov(1,1) - A_cov(0,1)*A_cov(1,0);

            libMesh::Real I3 = lambda_sq*det_A/det_a;

            libMesh::TensorValue<libMesh::Real> tau;
            this->_stress_strain_law.compute_stress(2,a_contra,a_cov,A_contra,A_cov,tau);

            if( _stress_indices[0] == quantity_index )
              {
                // Need to convert to Cauchy stress
                value = tau(0,0)/std::sqrt(I3);
              }
            else if( _stress_indices[1] == quantity_index )
              {
                // Need to convert to Cauchy stress
                value = tau(0,1)/std::sqrt(I3);
              }
            else if( _stress_indices[2] == quantity_index )
              {
                // Need to convert to Cauchy stress
                value = tau(1,1)/std::sqrt(I3);
              }
            else if( _stress_indices[3] == quantity_index )
              {
                value = 0.5*(tau(0,0) + tau(1,1)) + std::sqrt(0.25*(tau(0,0)-tau(1,1))*(tau(0,0)-tau(1,1))
                                                              + tau(0,1)*tau(0,1) );
              }
            else if( _stress_indices[4] == quantity_index ||
                     _stress_indices[5] == quantity_index ||
                     _stress_indices[6] == quantity_index   )
              {
                if(this->_is_compressible)
                  {
                    tau(2,2) = this->_stress_strain_law.compute_33_stress( a_contra, a_cov, A_contra, A_cov );
                  }

                libMesh::Real stress_I1 = tau(0,0) + tau(1,1) + tau(2,2);
                libMesh::Real stress_I2 = 0.5*(stress_I1*stress_I1 - (tau(0,0)*tau(0,0) + tau(1,1)*tau(1,1)
                                                                      + tau(2,2)*tau(2,2) + tau(0,1)*tau(0,1)
                                                                      + tau(1,0)*tau(1,0)) );

                libMesh::Real stress_I3 = tau(2,2)*( tau(0,0)*tau(1,1) - tau(1,0)*tau(0,1) );

                /* Formulae for principal stresses from:
                   http://en.wikiversity.org/wiki/Principal_stresses */

                // I_2^2 - 3*I_2
                libMesh::Real C1 = (stress_I1*stress_I1 - 3*stress_I2);

                // 2*I_1^3 - 9*I_1*_I2 + 27*I_3
                libMesh::Real C2 = (2*stress_I1*stress_I1*stress_I1 - 9*stress_I1*stress_I2 + 27*stress_I3)/54;

                libMesh::Real theta = std::acos( C2/(2*std::sqrt(C1*C1*C1)) )/3.0;

                if( _stress_indices[4] == quantity_index )
                  {
                    value = (stress_I1 + 2.0*std::sqrt(C1)*std::cos(theta))/3.0;
                  }

                if( _stress_indices[5] == quantity_index )
                  {
                    value = (stress_I1 + 2.0*std::sqrt(C1)*std::cos(theta+Constants::two_pi/3.0))/3.0;
                  }

                if( _stress_indices[6] == quantity_index )
                  {
                    value = (stress_I1 + 2.0*std::sqrt(C1)*std::cos(theta+2.0*Constants::two_pi/3.0))/3.0;
                  }
              }
            else if( _stress_zz_index == quantity_index )
              {
                value = this->_stress_strain_law.compute_33_stress( a_contra, a_cov, A_contra, A_cov );
              }
            else
              {
                //Wat?!
                libmesh_error();
              }
          } // is_stress

      }
  }

  template<typename StressStrainLaw>
  void ElasticMembrane<StressStrainLaw>::get_grad_disp( const AssemblyContext & context,
                                                        unsigned int qp,
                                                        libMesh::Gradient & grad_u,
                                                        libMesh::Gradient & grad_v,
                                                        libMesh::Gradient & grad_w )
  {
    const int n_u_dofs = context.get_dof_indices(this->_disp_vars.u()).size();

    // All shape function gradients are w.r.t. master element coordinates
    const std::vector<std::vector<libMesh::Real> > & dphi_dxi  = this->get_fe(context)->get_dphidxi();
    const std::vector<std::vector<libMesh::Real> > & dphi_deta = this->get_fe(context)->get_dphideta();

    const libMesh::DenseSubVector<libMesh::Number>& u_coeffs = context.get_elem_solution( this->_disp_vars.u() );
    const libMesh::DenseSubVector<libMesh::Number>& v_coeffs = context.get_elem_solution( this->_disp_vars.v() );
    const libMesh::DenseSubVector<libMesh::Number>* w_coeffs = NULL;

    if( this->_disp_vars.dim() == 3 )
      w_coeffs = &context.get_elem_solution( this->_disp_vars.w() );

    // Compute gradients  w.r.t. master element coordinates
    for( int d = 0; d < n_u_dofs; d++ )
      {
        libMesh::RealGradient u_gradphi( dphi_dxi[d][qp], dphi_deta[d][qp] );
        grad_u += u_coeffs(d)*u_gradphi;
        grad_v += v_coeffs(d)*u_gradphi;

        if( this->_disp_vars.dim() == 3 )
          grad_w += (*w_coeffs)(d)*u_gradphi;
      }
  }

  template<typename StressStrainLaw>
  void ElasticMembrane<StressStrainLaw>::get_stress_and_elasticity( const AssemblyContext & context,
                                                                    unsigned int qp,
                                                                    const libMesh::Gradient & grad_u,
                                                                    const libMesh::Gradient & grad_v,
                                                                    const libMesh::Gradient & grad_w,
                                                                    libMesh::TensorValue<libMesh::Real> & tau,
                                                                    ElasticityTensor & C )
  {
    libMesh::TensorValue<libMesh::Real> a_cov, a_contra, A_cov, A_contra;
    libMesh::Real lambda_sq = 0;

    this->compute_metric_tensors( qp, *(this->get_fe(context)), context,
                                  grad_u, grad_v, grad_w,
                                  a_cov, a_contra, A_cov, A_contra,
                                  lambda_sq );

    // Compute stress tensor
    const unsigned int dim = 2; // The membrane dimension is always 2 for this physics
    this->_stress_strain_law.compute_stress_and_elasticity(dim,a_contra,a_cov,A_contra,A_cov,tau,C);
  }

} // end namespace GRINS
