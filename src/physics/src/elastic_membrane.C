//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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
#include "grins/assembly_context.h"
#include "grins/solid_mechanics_bc_handling.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/quadrature.h"
#include "libmesh/boundary_info.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  template<typename StressStrainLaw>
  ElasticMembrane<StressStrainLaw>::ElasticMembrane( const GRINS::PhysicsName& physics_name, const GetPot& input,
                                                     bool lambda_sq_coupled, bool lambda_sq_var )
    : ElasticMembraneBase(physics_name,input),
      _stress_strain_law(input),
      _lambda_sq_coupled(lambda_sq_coupled),
      _lambda_sq_var(lambda_sq_var)
  {
    this->_bc_handler = new SolidMechanicsBCHandling( physics_name, input );

    return;
  }
  
  template<typename StressStrainLaw>
  ElasticMembrane<StressStrainLaw>::~ElasticMembrane()
  {
    return;
  }

  template<typename StressStrainLaw>
  void ElasticMembrane<StressStrainLaw>::element_time_derivative( bool compute_jacobian,
                                                                   AssemblyContext& context,
                                                                   CachedValues& /*cache*/ )
  {
    const unsigned int n_u_dofs = context.get_dof_indices(_disp_vars.u_var()).size();
    
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(_disp_vars.u_var())->get_JxW();

    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(_disp_vars.u_var())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(_disp_vars.u_var());
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(_disp_vars.v_var());
    libMesh::DenseSubVector<libMesh::Number> &Fw = context.get_elem_residual(_disp_vars.w_var());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    const std::vector<libMesh::RealGradient>& dxdxi  = context.get_element_fe(_disp_vars.u_var())->get_dxyzdxi();
    const std::vector<libMesh::RealGradient>& dxdeta = context.get_element_fe(_disp_vars.u_var())->get_dxyzdeta();

    const std::vector<libMesh::Real>& dxidx  = context.get_element_fe(_disp_vars.u_var())->get_dxidx();
    const std::vector<libMesh::Real>& dxidy  = context.get_element_fe(_disp_vars.u_var())->get_dxidy();
    const std::vector<libMesh::Real>& dxidz  = context.get_element_fe(_disp_vars.u_var())->get_dxidz();

    const std::vector<libMesh::Real>& detadx  = context.get_element_fe(_disp_vars.u_var())->get_detadx();
    const std::vector<libMesh::Real>& detady  = context.get_element_fe(_disp_vars.u_var())->get_detady();
    const std::vector<libMesh::Real>& detadz  = context.get_element_fe(_disp_vars.u_var())->get_detadz();

    const unsigned int dim = 2; // The manifold dimension is always 2 for this physics

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Gradient grad_u = context.interior_gradient(_disp_vars.u_var(), qp);
        libMesh::Gradient grad_v = context.interior_gradient(_disp_vars.v_var(), qp);
        libMesh::Gradient grad_w = context.interior_gradient(_disp_vars.w_var(), qp);

        libMesh::RealGradient grad_x( dxdxi[qp](0), dxdeta[qp](0) );
        libMesh::RealGradient grad_y( dxdxi[qp](1), dxdeta[qp](1) );
        libMesh::RealGradient grad_z( dxdxi[qp](2), dxdeta[qp](2) );

        libMesh::RealGradient dudxi( grad_u(0), grad_v(0), grad_w(0) );
        libMesh::RealGradient dudeta( grad_u(1), grad_v(1), grad_w(1) );
        
        libMesh::RealGradient dxi( dxidx[qp], dxidy[qp], dxidz[qp] );
        libMesh::RealGradient deta( detadx[qp], detady[qp], detadz[qp] );

        // Covariant metric tensor of reference configuration
        libMesh::TensorValue<libMesh::Real> a_cov( dxdxi[qp]*dxdxi[qp], dxdxi[qp]*dxdeta[qp], 0.0,
                                                   dxdeta[qp]*dxdxi[qp], dxdeta[qp]*dxdeta[qp] );

        // Covariant metric tensor of current configuration
        libMesh::TensorValue<libMesh::Real> A_cov( (dxdxi[qp] + dudxi)*(dxdxi[qp] + dudxi),
                                                   (dxdxi[qp] + dudxi)*(dxdeta[qp] + dudeta), 0.0,
                                                   (dxdeta[qp] + dudeta)*(dxdxi[qp] + dudxi),
                                                   (dxdeta[qp] + dudeta)*(dxdeta[qp] + dudeta) );

        // Contravariant metric tensor of reference configuration
        libMesh::TensorValue<libMesh::Real> a_contra( dxi*dxi, dxi*deta, 0.0,
                                                      deta*dxi, deta*deta );

        // Contravariant metric tensor in current configuration is A_cov^{-1}
        libMesh::Real det_A = A_cov(0,0)*A_cov(1,1) - A_cov(0,1)*A_cov(1,0);
        libMesh::TensorValue<libMesh::Real> A_contra(  A_cov(1,1)/det_A, -A_cov(0,1)/det_A, 0.0,
                                                      -A_cov(1,0)/det_A,  A_cov(0,0)/det_A );

        if( _lambda_sq_coupled )
          {
            a_cov(2,2)    = 1.0;
            a_contra(2,2) = 1.0;

            libMesh::Real det_a = a_cov.det();

            // If the material is incompressible, lambda^2 is known
            libMesh::Real lambda_sq = det_a/det_A;

            // If the material is compressible, then lambda_sq is an independent variable
            if( _lambda_sq_var )
              {
                libmesh_not_implemented();
              }

            A_cov(2,2) = lambda_sq;
            A_contra(2,2) = 1.0/lambda_sq;
          }

        // Compute stress tensor
        libMesh::TensorValue<libMesh::Real> tau;
        _stress_strain_law.compute_stress(dim,a_contra,a_cov,A_contra,A_cov,tau);

        libMesh::Real jac = JxW[qp];

        for (unsigned int i=0; i != n_u_dofs; i++)
	  {
            for( unsigned int alpha = 0; alpha < dim; alpha++ )
              {
                for( unsigned int beta = 0; beta < dim; beta++ )
                  {
                    Fu(i) -= 0.5*tau(alpha,beta)*( (grad_x(beta) + grad_u(beta))*u_gradphi[i][qp](alpha) +
                                                   (grad_x(alpha) + grad_u(alpha))*u_gradphi[i][qp](beta) )*jac;

                    Fv(i) -= 0.5*tau(alpha,beta)*( (grad_y(beta) + grad_v(beta))*u_gradphi[i][qp](alpha) +
                                                   (grad_y(alpha) + grad_v(alpha))*u_gradphi[i][qp](beta) )*jac;

                    Fw(i) -= 0.5*tau(alpha,beta)*( (grad_z(beta) + grad_w(beta))*u_gradphi[i][qp](alpha) +
                                                   (grad_z(alpha) + grad_w(alpha))*u_gradphi[i][qp](beta) )*jac;
                  }
              }
            if( compute_jacobian )
              {
                libmesh_not_implemented();
              }
          }

      }

    return;
  }

  template<typename StressStrainLaw>
  void ElasticMembrane<StressStrainLaw>::side_time_derivative( bool /*compute_jacobian*/,
                                                               AssemblyContext& /*context*/,
                                                               CachedValues& /*cache*/ )
  {
    /*
      std::vector<BoundaryID> ids = context.side_boundary_ids();
    
      for( std::vector<BoundaryID>::const_iterator it = ids.begin();
      it != ids.end(); it++ )
      {
      libmesh_assert (*it != libMesh::BoundaryInfo::invalid_id);
        
      _bc_handler->apply_neumann_bcs( context, cache, compute_jacobian, *it );
      } 
    */

    return;
  }

  template<typename StressStrainLaw>
  void ElasticMembrane<StressStrainLaw>::mass_residual( bool /*compute_jacobian*/,
                                                        AssemblyContext& /*context*/,
                                                        CachedValues& /*cache*/ )
  {
    libmesh_not_implemented();
    return;
  }

} // end namespace GRINS
