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
#include "grins/generic_ic_handler.h"
#include "grins/elasticity_tensor.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/quadrature.h"
#include "libmesh/boundary_info.h"
#include "libmesh/fem_system.h"
#include "libmesh/fe_interface.h"

namespace GRINS
{
  template<typename StressStrainLaw>
  ElasticMembrane<StressStrainLaw>::ElasticMembrane( const GRINS::PhysicsName& physics_name, const GetPot& input,
                                                     bool lambda_sq_coupled, bool lambda_sq_var )
    : ElasticMembraneBase(physics_name,input),
      _stress_strain_law(input),
      _h0( input("Physics/"+physics_name+"/h0", 1.0 ) ),
      _lambda_sq_coupled(lambda_sq_coupled),
      _lambda_sq_var(lambda_sq_var)
  {
    // Force the user to set h0
    if( !input.have_variable("Physics/"+physics_name+"/h0") )
      {
        std::cerr << "Error: Must specify initial thickness for "+physics_name << std::endl
                  << "       Input the option Physics/"+physics_name+"/h0" << std::endl;
        libmesh_error();
      }

    this->_bc_handler = new SolidMechanicsBCHandling( physics_name, input );

    this->_ic_handler = new GenericICHandler(physics_name, input);

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

    // Residuals that we're populating
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(_disp_vars.u_var());
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(_disp_vars.v_var());
    libMesh::DenseSubVector<libMesh::Number> &Fw = context.get_elem_residual(_disp_vars.w_var());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    // All shape function gradients are w.r.t. master element coordinates
    const std::vector<std::vector<libMesh::Real> >& dphi_dxi =
      context.get_element_fe(_disp_vars.u_var())->get_dphidxi();

    const std::vector<std::vector<libMesh::Real> >& dphi_deta =
      context.get_element_fe(_disp_vars.u_var())->get_dphideta();

    const libMesh::DenseSubVector<libMesh::Number>& u_coeffs = context.get_elem_solution( _disp_vars.u_var() );
    const libMesh::DenseSubVector<libMesh::Number>& v_coeffs = context.get_elem_solution( _disp_vars.v_var() );
    const libMesh::DenseSubVector<libMesh::Number>& w_coeffs = context.get_elem_solution( _disp_vars.w_var() );

    // Need these to build up the covariant and contravariant metric tensors
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

        libMesh::RealGradient dudxi( grad_u(0), grad_v(0), grad_w(0) );
        libMesh::RealGradient dudeta( grad_u(1), grad_v(1), grad_w(1) );
        
        libMesh::RealGradient dxi( dxidx[qp], dxidy[qp], dxidz[qp] );
        libMesh::RealGradient deta( detadx[qp], detady[qp], detadz[qp] );

        // Covariant metric tensor of reference configuration
        libMesh::TensorValue<libMesh::Real> a_cov( dxdxi[qp]*dxdxi[qp], dxdxi[qp]*dxdeta[qp], 0.0,
                                                   dxdeta[qp]*dxdxi[qp], dxdeta[qp]*dxdeta[qp] );

        libMesh::Real det_a = a_cov(0,0)*a_cov(1,1) - a_cov(0,1)*a_cov(1,0);

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
        //ElasticityTensor C;
        //_stress_strain_law.compute_stress_and_elasticity(dim,a_contra,a_cov,A_contra,A_cov,tau,C);

        libMesh::Real jac = JxW[qp];

        for (unsigned int i=0; i != n_u_dofs; i++)
	  {
            libMesh::RealGradient u_gradphi( dphi_dxi[i][qp], dphi_deta[i][qp] );

            for( unsigned int alpha = 0; alpha < dim; alpha++ )
              {
                for( unsigned int beta = 0; beta < dim; beta++ )
                  {
                    Fu(i) -= 0.5*tau(alpha,beta)*_h0*( (grad_x(beta) + grad_u(beta))*u_gradphi(alpha) +
                                                       (grad_x(alpha) + grad_u(alpha))*u_gradphi(beta) )*jac;

                    Fv(i) -= 0.5*tau(alpha,beta)*_h0*( (grad_y(beta) + grad_v(beta))*u_gradphi(alpha) +
                                                       (grad_y(alpha) + grad_v(alpha))*u_gradphi(beta) )*jac;

                    Fw(i) -= 0.5*tau(alpha,beta)*_h0*( (grad_z(beta) + grad_w(beta))*u_gradphi(alpha) +
                                                       (grad_z(alpha) + grad_w(alpha))*u_gradphi(beta) )*jac;
                  }
              }
          }

        if( compute_jacobian )
          {
            libmesh_not_implemented();
            /*
            for (unsigned int i=0; i != n_u_dofs; i++)
              {
                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    for( unsigned int alpha = 0; alpha < dim; alpha++ )
                      {
                        for( unsigned int beta = 0; beta < dim; beta++ )
                          {
                            for( unsigned int alpha = 0; alpha < dim; alpha++ )
                              {
                                for( unsigned int beta = 0; beta < dim; beta++ )
                                  {
                                    Kuu(i,j) -= 0.5*_h0*jac*( tau(alpha,beta)*( u_gradphi[j][qp](beta)*u_gradphi[i][qp](alpha) + u_gradphi[j][qp](alpha)*u_gradphi[i][qp](beta) )
                                                              + C(alpha,beta,lambda,mu)* *( (grad_x(beta) + grad_u(beta))*u_gradphi[i][qp](alpha) +
                                                                                            (grad_x(alpha) + grad_u(alpha))*u_gradphi[i][qp](beta) ) );
                                    Kuv(i,j) -= ;
                                    Kuw(i,j) -= ;

                                    Kvu(i,j) -= ;
                                    Kvv(i,j) -= ;
                                    Kvw(i,j) -= ;

                                    Kwu(i,j) -= ;
                                    Kwv(i,j) -= ;
                                    Kww(i,j) -= ;
                                  }
                              }
                          }
                      }
                  }
              }
            */
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

  template<typename StressStrainLaw>
  void ElasticMembrane<StressStrainLaw>::compute_postprocessed_quantity( unsigned int quantity_index,
                                                                         const AssemblyContext& context,
                                                                         const libMesh::Point& point,
                                                                         libMesh::Real& value )
  {
    value = std::numeric_limits<libMesh::Real>::quiet_NaN();

    bool is_stress = ( _stress_indices[0] == quantity_index ||
                       _stress_indices[1] == quantity_index ||
                       _stress_indices[2] == quantity_index   );

    bool is_strain = ( _strain_indices[0] == quantity_index ||
                       _strain_indices[1] == quantity_index ||
                       _strain_indices[2] == quantity_index   );

    if( is_stress || is_strain )
      {
        const unsigned int n_u_dofs = context.get_dof_indices(_disp_vars.u_var()).size();

        const libMesh::DenseSubVector<libMesh::Number>& u_coeffs = context.get_elem_solution( _disp_vars.u_var() );
        const libMesh::DenseSubVector<libMesh::Number>& v_coeffs = context.get_elem_solution( _disp_vars.v_var() );
        const libMesh::DenseSubVector<libMesh::Number>& w_coeffs = context.get_elem_solution( _disp_vars.w_var() );

        // Build new FE for the current point. We need this to build tensors at point.
        libMesh::AutoPtr<libMesh::FEGenericBase<libMesh::Real> > fe_new =
          this->build_new_fe( context.get_elem(), context.get_element_fe(_disp_vars.u_var()),
                              point );

        const std::vector<std::vector<libMesh::Real> >& dphi_dxi =
          fe_new->get_dphidxi();

        const std::vector<std::vector<libMesh::Real> >& dphi_deta =
          fe_new->get_dphideta();

        // Need these to build up the covariant and contravariant metric tensors
        const std::vector<libMesh::RealGradient>& dxdxi  = fe_new->get_dxyzdxi();
        const std::vector<libMesh::RealGradient>& dxdeta = fe_new->get_dxyzdeta();

        const std::vector<libMesh::Real>& dxidx  = fe_new->get_dxidx();
        const std::vector<libMesh::Real>& dxidy  = fe_new->get_dxidy();
        const std::vector<libMesh::Real>& dxidz  = fe_new->get_dxidz();

        const std::vector<libMesh::Real>& detadx  = fe_new->get_detadx();
        const std::vector<libMesh::Real>& detady  = fe_new->get_detady();
        const std::vector<libMesh::Real>& detadz  = fe_new->get_detadz();

        // Gradients are w.r.t. master element coordinates
        libMesh::Gradient grad_u, grad_v, grad_w;
        for( unsigned int d = 0; d < n_u_dofs; d++ )
          {
            libMesh::RealGradient u_gradphi( dphi_dxi[d][0], dphi_deta[d][0] );
            grad_u += u_coeffs(d)*u_gradphi;
            grad_v += v_coeffs(d)*u_gradphi;
            grad_w += w_coeffs(d)*u_gradphi;
          }

        libMesh::RealGradient grad_x( dxdxi[0](0), dxdeta[0](0) );
        libMesh::RealGradient grad_y( dxdxi[0](1), dxdeta[0](1) );
        libMesh::RealGradient grad_z( dxdxi[0](2), dxdeta[0](2) );

        libMesh::RealGradient dudxi( grad_u(0), grad_v(0), grad_w(0) );
        libMesh::RealGradient dudeta( grad_u(1), grad_v(1), grad_w(1) );

        libMesh::RealGradient dxi( dxidx[0], dxidy[0], dxidz[0] );
        libMesh::RealGradient deta( detadx[0], detady[0], detadz[0] );

        // Covariant metric tensor of reference configuration
        libMesh::TensorValue<libMesh::Real> a_cov( dxdxi[0]*dxdxi[0], dxdxi[0]*dxdeta[0], 0.0,
                                                   dxdeta[0]*dxdxi[0], dxdeta[0]*dxdeta[0] );

        libMesh::Real det_a = a_cov(0,0)*a_cov(1,1) - a_cov(0,1)*a_cov(1,0);

        // Covariant metric tensor of current configuration
        libMesh::TensorValue<libMesh::Real> A_cov( (dxdxi[0] + dudxi)*(dxdxi[0] + dudxi),
                                                   (dxdxi[0] + dudxi)*(dxdeta[0] + dudeta), 0.0,
                                                   (dxdeta[0] + dudeta)*(dxdxi[0] + dudxi),
                                                   (dxdeta[0] + dudeta)*(dxdeta[0] + dudeta) );

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

        // Contravariant metric tensor of reference configuration
        libMesh::TensorValue<libMesh::Real> a_contra( dxi*dxi, dxi*deta, 0.0,
                                                      deta*dxi, deta*deta );

        // Contravariant metric tensor in current configuration is A_cov^{-1}
        libMesh::Real det_A = A_cov(0,0)*A_cov(1,1) - A_cov(0,1)*A_cov(1,0);
        libMesh::TensorValue<libMesh::Real> A_contra(  A_cov(1,1)/det_A, -A_cov(0,1)/det_A, 0.0,
                                                      -A_cov(1,0)/det_A,  A_cov(0,0)/det_A );

        libMesh::Real I3 = det_A/det_a;

        if( _lambda_sq_coupled )
          {
            a_cov(2,2)    = 1.0;
            a_contra(2,2) = 1.0;

            // If the material is incompressible, lambda^2 is known
            libMesh::Real lambda_sq = det_a/det_A;

            // If the material is compressible, then lambda_sq is an independent variable
            if( _lambda_sq_var )
              {
                libmesh_not_implemented();
              }

            A_cov(2,2) = lambda_sq;
            A_contra(2,2) = 1.0/lambda_sq;

            I3 *= lambda_sq;
          }

        libMesh::TensorValue<libMesh::Real> tau;
        _stress_strain_law.compute_stress(2,a_contra,a_cov,A_contra,A_cov,tau);

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
        else
          {
            //Wat?!
            libmesh_error();
          }

      }

    return;
  }

  template<typename StressStrainLaw>
  libMesh::AutoPtr<libMesh::FEGenericBase<libMesh::Real> > ElasticMembrane<StressStrainLaw>::build_new_fe( const libMesh::Elem& elem, 
                                                                                                           const libMesh::FEGenericBase<libMesh::Real>* fe,
                                                                                                           const libMesh::Point p )
  {
    using namespace libMesh;
    FEType fe_type = fe->get_fe_type();

    // If we don't have an Elem to evaluate on, then the only functions
    // we can sensibly evaluate are the scalar dofs which are the same
    // everywhere.
    libmesh_assert(&elem || fe_type.family == SCALAR);

    unsigned int elem_dim = &elem ? elem.dim() : 0;

    AutoPtr<FEGenericBase<libMesh::Real> >
      fe_new(FEGenericBase<libMesh::Real>::build(elem_dim, fe_type));

    // Map the physical co-ordinates to the master co-ordinates using the inverse_map from fe_interface.h
    // Build a vector of point co-ordinates to send to reinit
    Point master_point = &elem ?
      FEInterface::inverse_map(elem_dim, fe_type, &elem, p) :
      Point(0);

    std::vector<Point> coor(1, master_point);

    // Reinitialize the element and compute the shape function values at coor
    fe_new->reinit (&elem, &coor);

    return fe_new;
  }

} // end namespace GRINS
