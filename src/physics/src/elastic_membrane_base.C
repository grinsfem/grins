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
#include "grins/elastic_membrane_base.h"

// GRINS
#include "grins/materials_parsing.h"
#include "grins/multiphysics_sys.h"

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
                                     "MembraneThickness",
                                     PhysicsNaming::elastic_membrane(),
                                     (*this),
                                     _h0 );

    MaterialsParsing::read_property( input,
                                     "Density",
                                     PhysicsNaming::elastic_membrane(),
                                     (*this),
                                     _rho );

    if( this->_disp_vars.dim() < 2 )
      libmesh_error_msg("ERROR: ElasticMembraneBase subclasses only valid for two or three dimensions! Make sure you have at least two components in your Displacement type variable.");
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

    const MultiphysicsSystem & system = context.get_multiphysics_system();

    unsigned int u_dot_var = system.get_second_order_dot_var(this->_disp_vars.u());
    unsigned int v_dot_var = system.get_second_order_dot_var(this->_disp_vars.v());

    // Residuals that we're populating
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(u_dot_var);
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(v_dot_var);
    libMesh::DenseSubVector<libMesh::Number>* Fw = NULL;

    libMesh::DenseSubMatrix<libMesh::Number>& Kuu = context.get_elem_jacobian(u_dot_var,u_dot_var);
    libMesh::DenseSubMatrix<libMesh::Number>& Kvv = context.get_elem_jacobian(v_dot_var,v_dot_var);
    libMesh::DenseSubMatrix<libMesh::Number>* Kww = NULL;

    unsigned int w_dot_var = libMesh::invalid_uint;
    if( this->_disp_vars.dim() == 3 )
      {
        w_dot_var = system.get_second_order_dot_var(this->_disp_vars.w());
        Fw = &context.get_elem_residual(w_dot_var);
        Kww = &context.get_elem_jacobian(w_dot_var,w_dot_var);
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Real jac = JxW[qp];

        libMesh::Real u_ddot, v_ddot;
        (context.*interior_solution)( u_dot_var, qp, u_ddot );
        (context.*interior_solution)( v_dot_var, qp, v_ddot );

        libMesh::Real w_ddot = 0.0;
        if( this->_disp_vars.dim() == 3 )
          (context.*interior_solution)( w_dot_var, qp, w_ddot );

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) += mu*this->_rho*_h0*u_ddot*u_phi[i][qp]*jac;
            Fv(i) += mu*this->_rho*_h0*v_ddot*u_phi[i][qp]*jac;

            if( this->_disp_vars.dim() == 3 )
              (*Fw)(i) += mu*this->_rho*_h0*w_ddot*u_phi[i][qp]*jac;

            if( compute_jacobian )
              {
                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    libMesh::Real jac_term = this->_rho*_h0*u_phi[i][qp]*u_phi[j][qp]*jac;
                    jac_term *= mu*(context.*get_solution_deriv)();

                    Kuu(i,j) += jac_term;
                    Kvv(i,j) += jac_term;

                    if( this->_disp_vars.dim() == 3 )
                      (*Kww)(i,j) += jac_term;
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
