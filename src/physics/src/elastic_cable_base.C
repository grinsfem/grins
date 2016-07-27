//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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
  void ElasticCableBase<StressStrainLaw>::mass_residual_impl( bool compute_jacobian,
                                                              AssemblyContext& context,
                                                              InteriorFuncType interior_solution,
                                                              VarDerivType get_solution_deriv,
                                                              libMesh::Real mu )
  {
    const unsigned int n_u_dofs = context.get_dof_indices(this->_disp_vars.u()).size();

    const std::vector<libMesh::Real> &JxW =
      this->get_fe(context)->get_JxW();

    const std::vector<std::vector<libMesh::Real> >& u_phi =
      this->get_fe(context)->get_phi();

    // Residuals that we're populating
    libMesh::DenseSubVector<libMesh::Number> & Fu = context.get_elem_residual(this->_disp_vars.u());
    libMesh::DenseSubVector<libMesh::Number> * Fv = NULL;
    libMesh::DenseSubVector<libMesh::Number> * Fw = NULL;

    libMesh::DenseSubMatrix<libMesh::Number> & Kuu = context.get_elem_jacobian(this->_disp_vars.u(),this->_disp_vars.u());
    libMesh::DenseSubMatrix<libMesh::Number> * Kvv = NULL;
    libMesh::DenseSubMatrix<libMesh::Number> * Kww = NULL;

    if( this->_disp_vars.dim() >= 2 )
      {
        Fv = &context.get_elem_residual(this->_disp_vars.v());
        Kvv = &context.get_elem_jacobian(this->_disp_vars.v(),this->_disp_vars.v());
      }

    if( this->_disp_vars.dim() == 3 )
      {
        Fw = &context.get_elem_residual(this->_disp_vars.w());
        Kww = &context.get_elem_jacobian(this->_disp_vars.w(),this->_disp_vars.w());
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Real jac = JxW[qp];

        libMesh::Real u_ddot, v_ddot, w_ddot;
        (context.*interior_solution)( this->_disp_vars.u(), qp, u_ddot );

        if( this->_disp_vars.dim() >= 2 )
          (context.*interior_solution)( this->_disp_vars.v(), qp, v_ddot );

        if( this->_disp_vars.dim() == 3 )
          (context.*interior_solution)( this->_disp_vars.w(), qp, w_ddot );

        for (unsigned int i=0; i != n_u_dofs; i++)
	  {
            libMesh::Real value = this->_rho*this->_A*u_phi[i][qp]*jac*mu;
            Fu(i) += value*u_ddot;

            if( this->_disp_vars.dim() >= 2 )
              (*Fv)(i) += value*v_ddot;

            if( this->_disp_vars.dim() == 3 )
              (*Fw)(i) += value*w_ddot;

            if( compute_jacobian )
              {
                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    libMesh::Real jac_term = mu*this->_rho*this->_A*u_phi[i][qp]*u_phi[j][qp]*jac;
                    jac_term *= (context.*get_solution_deriv)();

                    Kuu(i,j) += jac_term;

                    if( this->_disp_vars.dim() >= 2 )
                      (*Kvv)(i,j) += jac_term;

                    if( this->_disp_vars.dim() == 3 )
                      (*Kww)(i,j) += jac_term;
                  }
              }
          }
      }
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
