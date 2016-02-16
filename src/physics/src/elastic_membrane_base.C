//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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
#include "grins/elastic_membrane_base.h"

// GRINS
#include "grins/elasticity_tensor.h"
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/quadrature.h"
#include "libmesh/fem_system.h"
#include "libmesh/elem.h"

namespace GRINS
{
  template<typename StressStrainLaw>
  ElasticMembraneBase<StressStrainLaw>::ElasticMembraneBase( const GRINS::PhysicsName& physics_name,
                                                             const GetPot& input,
                                                             bool is_compressible )
    : ElasticMembraneAbstract(physics_name,input),
      _stress_strain_law(input,MaterialsParsing::material_name(input,PhysicsNaming::elastic_membrane())),
      _is_compressible(is_compressible),
      _h0(0.0)
  {
    MaterialsParsing::read_property( input,
                                     "Physics/"+physics_name+"/h0",
                                     "MembraneThickness",
                                     PhysicsNaming::elastic_membrane(),
                                     (*this),
                                     _h0 );

    MaterialsParsing::read_density( physics_name,
                                    input,
                                    (*this),
                                    _rho );
  }

  template<typename StressStrainLaw>
  void ElasticMembraneBase<StressStrainLaw>::init_variables( libMesh::FEMSystem* system )
  {
    // First call base class
    ElasticMembraneAbstract::init_variables(system);

    // Now build lambda_sq variable if we need it
    if(_is_compressible)
      {
        /*! \todo Might want to make Order/FEType inputable */
        _lambda_sq_var = system->add_variable( "lambda_sq", GRINSEnums::FIRST, GRINSEnums::LAGRANGE);
      }
  }

  template<typename StressStrainLaw>
  void ElasticMembraneBase<StressStrainLaw>::element_time_derivative_impl( bool compute_jacobian,
                                                                           AssemblyContext& context,
                                                                           VarFuncType get_solution,
                                                                           VarDerivType get_solution_deriv,
                                                                           libMesh::Real lambda_factor )
  {
    const unsigned int n_u_dofs = context.get_dof_indices(_disp_vars.u()).size();

    const std::vector<libMesh::Real> &JxW =
      this->get_fe(context)->get_JxW();

    // Residuals that we're populating
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(_disp_vars.u());
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(_disp_vars.v());
    libMesh::DenseSubVector<libMesh::Number> &Fw = context.get_elem_residual(_disp_vars.w());

    libMesh::DenseSubMatrix<libMesh::Number>& Kuu = context.get_elem_jacobian(_disp_vars.u(),_disp_vars.u());
    libMesh::DenseSubMatrix<libMesh::Number>& Kuv = context.get_elem_jacobian(_disp_vars.u(),_disp_vars.v());
    libMesh::DenseSubMatrix<libMesh::Number>& Kuw = context.get_elem_jacobian(_disp_vars.u(),_disp_vars.w());

    libMesh::DenseSubMatrix<libMesh::Number>& Kvu = context.get_elem_jacobian(_disp_vars.v(),_disp_vars.u());
    libMesh::DenseSubMatrix<libMesh::Number>& Kvv = context.get_elem_jacobian(_disp_vars.v(),_disp_vars.v());
    libMesh::DenseSubMatrix<libMesh::Number>& Kvw = context.get_elem_jacobian(_disp_vars.v(),_disp_vars.w());

    libMesh::DenseSubMatrix<libMesh::Number>& Kwu = context.get_elem_jacobian(_disp_vars.w(),_disp_vars.u());
    libMesh::DenseSubMatrix<libMesh::Number>& Kwv = context.get_elem_jacobian(_disp_vars.w(),_disp_vars.v());
    libMesh::DenseSubMatrix<libMesh::Number>& Kww = context.get_elem_jacobian(_disp_vars.w(),_disp_vars.w());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    // All shape function gradients are w.r.t. master element coordinates
    const std::vector<std::vector<libMesh::Real> >& dphi_dxi =
      this->get_fe(context)->get_dphidxi();

    const std::vector<std::vector<libMesh::Real> >& dphi_deta =
      this->get_fe(context)->get_dphideta();

    const libMesh::DenseSubVector<libMesh::Number>& u_coeffs = (context.*get_solution)( _disp_vars.u() );
    const libMesh::DenseSubVector<libMesh::Number>& v_coeffs = (context.*get_solution)( _disp_vars.v() );
    const libMesh::DenseSubVector<libMesh::Number>& w_coeffs = (context.*get_solution)( _disp_vars.w() );

    // Need these to build up the covariant and contravariant metric tensors
    const std::vector<libMesh::RealGradient>& dxdxi  = this->get_fe(context)->get_dxyzdxi();
    const std::vector<libMesh::RealGradient>& dxdeta = this->get_fe(context)->get_dxyzdeta();

    const unsigned int dim = 2; // The manifold dimension is always 2 for this physics

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // Gradients are w.r.t. master element coordinates
        libMesh::Gradient grad_u, grad_v, grad_w;
        for( unsigned int d = 0; d < n_u_dofs; d++ )
          {
            libMesh::RealGradient u_gradphi( dphi_dxi[d][qp], dphi_deta[d][qp] );
            grad_u += u_coeffs(d)*u_gradphi;
            grad_v += v_coeffs(d)*u_gradphi;
            grad_w += w_coeffs(d)*u_gradphi;
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

        // Compute stress and elasticity tensors
        libMesh::TensorValue<libMesh::Real> tau;
        ElasticityTensor C;
        _stress_strain_law.compute_stress_and_elasticity(dim,a_contra,a_cov,A_contra,A_cov,tau,C);

        libMesh::Real jac = JxW[qp];

        for (unsigned int i=0; i != n_u_dofs; i++)
	  {
            libMesh::RealGradient u_gradphi( dphi_dxi[i][qp], dphi_deta[i][qp] );

            for( unsigned int alpha = 0; alpha < dim; alpha++ )
              {
                for( unsigned int beta = 0; beta < dim; beta++ )
                  {
                    libMesh::Real factor = lambda_factor*0.5*tau(alpha,beta)*_h0*jac;

                    Fu(i) += factor*( (grad_x(beta) + grad_u(beta))*u_gradphi(alpha) +
                                      (grad_x(alpha) + grad_u(alpha))*u_gradphi(beta) );

                    Fv(i) += factor*( (grad_y(beta) + grad_v(beta))*u_gradphi(alpha) +
                                      (grad_y(alpha) + grad_v(alpha))*u_gradphi(beta) );

                    Fw(i) += factor*( (grad_z(beta) + grad_w(beta))*u_gradphi(alpha) +
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

                    for( unsigned int alpha = 0; alpha < dim; alpha++ )
                      {
                        for( unsigned int beta = 0; beta < dim; beta++ )
                          {
                            const libMesh::Real diag_term = lambda_factor*0.5*_h0*jac*tau(alpha,beta)*(context.*get_solution_deriv)()*
                                                            ( u_gradphi_j(beta)*u_gradphi_i(alpha) +
                                                              u_gradphi_j(alpha)*u_gradphi_i(beta) );
                            Kuu(i,j) += diag_term;

                            Kvv(i,j) += diag_term;

                            Kww(i,j) += diag_term;

                            for( unsigned int lambda = 0; lambda < dim; lambda++ )
                              {
                                for( unsigned int mu = 0; mu < dim; mu++ )
                                  {
                                    const libMesh::Real dgamma_du = 0.5*( u_gradphi_j(lambda)*(grad_x(mu)+grad_u(mu)) +
                                                                          (grad_x(lambda)+grad_u(lambda))*u_gradphi_j(mu) );

                                    const libMesh::Real dgamma_dv = 0.5*( u_gradphi_j(lambda)*(grad_y(mu)+grad_v(mu)) +
                                                                          (grad_y(lambda)+grad_v(lambda))*u_gradphi_j(mu) );

                                    const libMesh::Real dgamma_dw = 0.5*( u_gradphi_j(lambda)*(grad_z(mu)+grad_w(mu)) +
                                                                          (grad_z(lambda)+grad_w(lambda))*u_gradphi_j(mu) );

                                    const libMesh::Real C1 = lambda_factor*0.5*_h0*jac*C(alpha,beta,lambda,mu)*(context.*get_solution_deriv)();

                                    const libMesh::Real x_term = C1*( (grad_x(beta)+grad_u(beta))*u_gradphi_i(alpha) +
                                                                      (grad_x(alpha)+grad_u(alpha))*u_gradphi_i(beta) );

                                    const libMesh::Real y_term = C1*( (grad_y(beta)+grad_v(beta))*u_gradphi_i(alpha) +
                                                                      (grad_y(alpha)+grad_v(alpha))*u_gradphi_i(beta) );

                                    const libMesh::Real z_term = C1*( (grad_z(beta)+grad_w(beta))*u_gradphi_i(alpha) +
                                                                      (grad_z(alpha)+grad_w(alpha))*u_gradphi_i(beta) );

                                    Kuu(i,j) += x_term*dgamma_du;

                                    Kuv(i,j) += x_term*dgamma_dv;

                                    Kuw(i,j) += x_term*dgamma_dw;

                                    Kvu(i,j) += y_term*dgamma_du;

                                    Kvv(i,j) += y_term*dgamma_dv;

                                    Kvw(i,j) += y_term*dgamma_dw;

                                    Kwu(i,j) += z_term*dgamma_du;

                                    Kwv(i,j) += z_term*dgamma_dv;

                                    Kww(i,j) += z_term*dgamma_dw;
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
  void ElasticMembraneBase<StressStrainLaw>::mass_residual_impl( bool compute_jacobian,
                                                                 AssemblyContext& context,
                                                                 InteriorFuncType interior_solution,
                                                                 VarDerivType get_solution_deriv,
                                                                 libMesh::Real mu )
  {
    const unsigned int n_u_dofs = context.get_dof_indices(_disp_vars.u()).size();

    const std::vector<libMesh::Real> &JxW =
      this->get_fe(context)->get_JxW();

    const std::vector<std::vector<libMesh::Real> >& u_phi =
      this->get_fe(context)->get_phi();

    // Residuals that we're populating
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(_disp_vars.u());
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(_disp_vars.v());
    libMesh::DenseSubVector<libMesh::Number> &Fw = context.get_elem_residual(_disp_vars.w());

    libMesh::DenseSubMatrix<libMesh::Number>& Kuu = context.get_elem_jacobian(_disp_vars.u(),_disp_vars.u());
    libMesh::DenseSubMatrix<libMesh::Number>& Kvv = context.get_elem_jacobian(_disp_vars.v(),_disp_vars.v());
    libMesh::DenseSubMatrix<libMesh::Number>& Kww = context.get_elem_jacobian(_disp_vars.w(),_disp_vars.w());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Real jac = JxW[qp];

        libMesh::Real u_ddot, v_ddot, w_ddot;
        (context.*interior_solution)( _disp_vars.u(), qp, u_ddot );
        (context.*interior_solution)( _disp_vars.v(), qp, v_ddot );
        (context.*interior_solution)( _disp_vars.w(), qp, w_ddot );

        for (unsigned int i=0; i != n_u_dofs; i++)
	  {
            Fu(i) += mu*this->_rho*_h0*u_ddot*u_phi[i][qp]*jac;
            Fv(i) += mu*this->_rho*_h0*v_ddot*u_phi[i][qp]*jac;
            Fw(i) += mu*this->_rho*_h0*w_ddot*u_phi[i][qp]*jac;

            if( compute_jacobian )
              {
                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    libMesh::Real jac_term = this->_rho*_h0*u_phi[i][qp]*u_phi[j][qp]*jac;
                    jac_term *= mu*(context.*get_solution_deriv)();

                    Kuu(i,j) += jac_term;
                    Kvv(i,j) += jac_term;
                    Kww(i,j) += jac_term;
                  }
              }
          }
      }
  }

  template<typename StressStrainLaw>
  void ElasticMembraneBase<StressStrainLaw>::compute_metric_tensors( unsigned int qp,
                                                                     const libMesh::FEBase& elem,
                                                                     const AssemblyContext& context,
                                                                     const libMesh::Gradient& grad_u,
                                                                     const libMesh::Gradient& grad_v,
                                                                     const libMesh::Gradient& grad_w,
                                                                     libMesh::TensorValue<libMesh::Real>& a_cov,
                                                                     libMesh::TensorValue<libMesh::Real>& a_contra,
                                                                     libMesh::TensorValue<libMesh::Real>& A_cov,
                                                                     libMesh::TensorValue<libMesh::Real>& A_contra,
                                                                     libMesh::Real& lambda_sq )
  {
    const std::vector<libMesh::RealGradient>& dxdxi  = elem.get_dxyzdxi();
    const std::vector<libMesh::RealGradient>& dxdeta = elem.get_dxyzdeta();

    const std::vector<libMesh::Real>& dxidx  = elem.get_dxidx();
    const std::vector<libMesh::Real>& dxidy  = elem.get_dxidy();
    const std::vector<libMesh::Real>& dxidz  = elem.get_dxidz();

    const std::vector<libMesh::Real>& detadx  = elem.get_detadx();
    const std::vector<libMesh::Real>& detady  = elem.get_detady();
    const std::vector<libMesh::Real>& detadz  = elem.get_detadz();

    libMesh::RealGradient dxi( dxidx[qp], dxidy[qp], dxidz[qp] );
    libMesh::RealGradient deta( detadx[qp], detady[qp], detadz[qp] );

    libMesh::RealGradient dudxi( grad_u(0), grad_v(0), grad_w(0) );
    libMesh::RealGradient dudeta( grad_u(1), grad_v(1), grad_w(1) );

    // Covariant metric tensor of reference configuration
    a_cov.zero();
    a_cov(0,0) = dxdxi[qp]*dxdxi[qp];
    a_cov(0,1) = dxdxi[qp]*dxdeta[qp];
    a_cov(1,0) = dxdeta[qp]*dxdxi[qp];
    a_cov(1,1) = dxdeta[qp]*dxdeta[qp];

    libMesh::Real det_a = a_cov(0,0)*a_cov(1,1) - a_cov(0,1)*a_cov(1,0);

    // Covariant metric tensor of current configuration
    A_cov.zero();
    A_cov(0,0) = (dxdxi[qp] + dudxi)*(dxdxi[qp] + dudxi);
    A_cov(0,1) = (dxdxi[qp] + dudxi)*(dxdeta[qp] + dudeta);
    A_cov(1,0) = (dxdeta[qp] + dudeta)*(dxdxi[qp] + dudxi);
    A_cov(1,1) = (dxdeta[qp] + dudeta)*(dxdeta[qp] + dudeta);

    // Contravariant metric tensor of reference configuration
    a_contra.zero();
    a_contra(0,0) = dxi*dxi;
    a_contra(0,1) = dxi*deta;
    a_contra(1,0) = deta*dxi;
    a_contra(1,1) = deta*deta;

    // Contravariant metric tensor in current configuration is A_cov^{-1}
    libMesh::Real det_A = A_cov(0,0)*A_cov(1,1) - A_cov(0,1)*A_cov(1,0);

    A_contra.zero();
    A_contra(0,0) =  A_cov(1,1)/det_A;
    A_contra(0,1) = -A_cov(0,1)/det_A;
    A_contra(1,0) = -A_cov(1,0)/det_A;
    A_contra(1,1) =  A_cov(0,0)/det_A;

    a_cov(2,2)    = 1.0;
    a_contra(2,2) = 1.0;


    // If the material is compressible, then lambda_sq is an independent variable
    if( _is_compressible )
      {
        lambda_sq = context.interior_value(this->_lambda_sq_var, qp);
      }
    else
      {
        // If the material is incompressible, lambda^2 is known
        lambda_sq = det_a/det_A;
      }

    A_cov(2,2) = lambda_sq;
    A_contra(2,2) = 1.0/lambda_sq;
  }

} // end namespace GRINS
