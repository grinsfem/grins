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
#include "grins/elastic_cable_base.h"

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
  ElasticCableBase<StressStrainLaw>::ElasticCableBase( const PhysicsName& physics_name,
                                                       const GetPot& input,
                                                       bool is_compressible)
    : ElasticCableAbstract(physics_name,input),
      _stress_strain_law(input,MaterialsParsing::material_name(input,PhysicsNaming::elastic_cable())),
      _is_compressible(is_compressible)
  {}

  template<typename StressStrainLaw>
  void ElasticCableBase<StressStrainLaw>::element_time_derivative_impl( bool compute_jacobian,
                                                                        AssemblyContext& context,
                                                                        VarFuncType get_solution,
                                                                        VarDerivType get_solution_deriv,
                                                                        libMesh::Real lambda )
  {
    const unsigned int n_u_dofs = context.get_dof_indices(_disp_vars.u()).size();

    const std::vector<libMesh::Real> &JxW =
      this->get_fe(context)->get_JxW();

    // Residuals that we're populating
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(_disp_vars.u());
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(_disp_vars.v());
    libMesh::DenseSubVector<libMesh::Number> &Fw = context.get_elem_residual(_disp_vars.w());

    //Grab the Jacobian matrix as submatrices
    //libMesh::DenseMatrix<libMesh::Number> &K = context.get_elem_jacobian();
    libMesh::DenseSubMatrix<libMesh::Number> &Kuu = context.get_elem_jacobian(_disp_vars.u(),_disp_vars.u());
    libMesh::DenseSubMatrix<libMesh::Number> &Kuv = context.get_elem_jacobian(_disp_vars.u(),_disp_vars.v());
    libMesh::DenseSubMatrix<libMesh::Number> &Kuw = context.get_elem_jacobian(_disp_vars.u(),_disp_vars.w());
    libMesh::DenseSubMatrix<libMesh::Number> &Kvu = context.get_elem_jacobian(_disp_vars.v(),_disp_vars.u());
    libMesh::DenseSubMatrix<libMesh::Number> &Kvv = context.get_elem_jacobian(_disp_vars.v(),_disp_vars.v());
    libMesh::DenseSubMatrix<libMesh::Number> &Kvw = context.get_elem_jacobian(_disp_vars.v(),_disp_vars.w());
    libMesh::DenseSubMatrix<libMesh::Number> &Kwu = context.get_elem_jacobian(_disp_vars.w(),_disp_vars.u());
    libMesh::DenseSubMatrix<libMesh::Number> &Kwv = context.get_elem_jacobian(_disp_vars.w(),_disp_vars.v());
    libMesh::DenseSubMatrix<libMesh::Number> &Kww = context.get_elem_jacobian(_disp_vars.w(),_disp_vars.w());


    unsigned int n_qpoints = context.get_element_qrule().n_points();

    // All shape function gradients are w.r.t. master element coordinates
    const std::vector<std::vector<libMesh::Real> >& dphi_dxi = this->get_fe(context)->get_dphidxi();

    const libMesh::DenseSubVector<libMesh::Number>& u_coeffs = (context.*get_solution)( _disp_vars.u() );
    const libMesh::DenseSubVector<libMesh::Number>& v_coeffs = (context.*get_solution)( _disp_vars.v() );
    const libMesh::DenseSubVector<libMesh::Number>& w_coeffs = (context.*get_solution)( _disp_vars.w() );

    // Need these to build up the covariant and contravariant metric tensors
    const std::vector<libMesh::RealGradient>& dxdxi  = this->get_fe(context)->get_dxyzdxi();

    const unsigned int dim = 1; // The cable dimension is always 1 for this physics

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // Gradients are w.r.t. master element coordinates
        libMesh::Gradient grad_u, grad_v, grad_w;

        for( unsigned int d = 0; d < n_u_dofs; d++ )
          {
            libMesh::RealGradient u_gradphi( dphi_dxi[d][qp] );
            grad_u += u_coeffs(d)*u_gradphi;
            grad_v += v_coeffs(d)*u_gradphi;
            grad_w += w_coeffs(d)*u_gradphi;
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
        //_stress_strain_law.compute_stress(dim,a_contra,a_cov,A_contra,A_cov,tau);
        _stress_strain_law.compute_stress_and_elasticity(dim,a_contra,a_cov,A_contra,A_cov,tau,C);


        libMesh::Real jac = JxW[qp];

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            libMesh::RealGradient u_gradphi( dphi_dxi[i][qp] );

            const libMesh::Real res_term = lambda*tau(0,0)*_A*jac*u_gradphi(0);

            Fu(i) += res_term*(grad_x(0) + grad_u(0));

            Fv(i) += res_term*(grad_y(0) + grad_v(0));

            Fw(i) += res_term*(grad_z(0) + grad_w(0));
          }

        if( compute_jacobian )
          {
            for(unsigned int i=0; i != n_u_dofs; i++)
              {
                libMesh::RealGradient u_gradphi_I( dphi_dxi[i][qp] );
                for(unsigned int j=0; j != n_u_dofs; j++)
                  {
                    libMesh::RealGradient u_gradphi_J( dphi_dxi[j][qp] );

                    const libMesh::Real diag_term = lambda*_A*jac*tau(0,0)*( u_gradphi_J(0)*u_gradphi_I(0))*(context.*get_solution_deriv)();

                    Kuu(i,j) += diag_term;

                    Kvv(i,j) += diag_term;

                    Kww(i,j) += diag_term;

                    const libMesh::Real dgamma_du = ( u_gradphi_J(0)*(grad_x(0)+grad_u(0)) );

                    const libMesh::Real dgamma_dv = ( u_gradphi_J(0)*(grad_y(0)+grad_v(0)) );

                    const libMesh::Real dgamma_dw = ( u_gradphi_J(0)*(grad_z(0)+grad_w(0)) );

                    const libMesh::Real C1 = lambda*_A*jac*C(0,0,0,0)*(context.*get_solution_deriv)();

                    const libMesh::Real x_term = C1*( (grad_x(0)+grad_u(0))*u_gradphi_I(0) );

                    const libMesh::Real y_term = C1*( (grad_y(0)+grad_v(0))*u_gradphi_I(0) );

                    const libMesh::Real z_term = C1*( (grad_z(0)+grad_w(0))*u_gradphi_I(0) );

                    Kuu(i,j) += x_term*dgamma_du;

                    Kuv(i,j) += x_term*dgamma_dv;

                    Kuw(i,j) += x_term*dgamma_dw;

                    Kvu(i,j) += y_term*dgamma_du;

                    Kvv(i,j) += y_term*dgamma_dv;

                    Kvw(i,j) += y_term*dgamma_dw;

                    Kwu(i,j) += z_term*dgamma_du;

                    Kwv(i,j) += z_term*dgamma_dv;

                    Kww(i,j) += z_term*dgamma_dw;

                  } // end j-loop
              } // end i-loop
          } // end if(compute_jacobian)
      } // end qp loop
  }

  template<typename StressStrainLaw>
  void ElasticCableBase<StressStrainLaw>::compute_metric_tensors( unsigned int qp,
                                                                  const libMesh::FEBase& elem,
                                                                  const AssemblyContext& /*context*/,
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

    const std::vector<libMesh::Real>& dxidx  = elem.get_dxidx();
    const std::vector<libMesh::Real>& dxidy  = elem.get_dxidy();
    const std::vector<libMesh::Real>& dxidz  = elem.get_dxidz();

    libMesh::RealGradient dxi( dxidx[qp], dxidy[qp], dxidz[qp] );

    libMesh::RealGradient dudxi( grad_u(0), grad_v(0), grad_w(0) );

    // Covariant metric tensor of reference configuration
    a_cov.zero();
    a_cov(0,0) = dxdxi[qp]*dxdxi[qp];
    a_cov(1,1)    = 1.0;
    a_cov(2,2)    = 1.0;

    // Covariant metric tensor of current configuration
    A_cov.zero();
    A_cov(0,0) = (dxdxi[qp] + dudxi)*(dxdxi[qp] + dudxi);

    // Contravariant metric tensor of reference configuration
    a_contra.zero();
    a_contra(0,0) = 1/a_cov(0,0);
    a_contra(1,1) = 1.0;
    a_contra(2,2) = 1.0;

    // Contravariant metric tensor in current configuration is A_cov^{-1}
    A_contra.zero();
    A_contra(0,0) =  1/A_cov(0,0);

    // If the material is compressible, then lambda_sq is an independent variable
    if( _is_compressible )
      {
        libmesh_not_implemented();
        //lambda_sq = context.interior_value(this->_lambda_sq_var, qp);
      }
    else
      {
        // If the material is incompressible, lambda^2 is known
        lambda_sq = a_cov(0,0)/A_cov(0,0);//det_a/det_A;
      }

    //Update the covariant and contravariant tensors of current configuration
    A_cov(1,1)    = lambda_sq;
    A_cov(2,2)    = lambda_sq;
    A_contra(1,1) = 1.0/lambda_sq;
    A_contra(2,2) = 1.0/lambda_sq;
  }

} // end namespace GRINS
