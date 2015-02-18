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
#include "grins/elastic_cable.h"

// GRINS
#include "grins_config.h"
#include "grins/assembly_context.h"
#include "grins/solid_mechanics_bc_handling.h"
#include "grins/generic_ic_handler.h"
#include "grins/elasticity_tensor.h"
#include "grins/postprocessed_quantities.h"
#include "grins/physics.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/quadrature.h"
#include "libmesh/boundary_info.h"
#include "libmesh/fem_system.h"
#include "libmesh/fe_interface.h"

namespace GRINS
{
  template<typename StressStrainLaw>
  ElasticCable<StressStrainLaw>::ElasticCable( const GRINS::PhysicsName& physics_name, const GetPot& input,
                                               bool is_compressible )
    : Physics(physics_name,input),
      _disp_vars(input,physics_name),
      _stress_strain_law(input),
      _A( input("Physics/"+physics_name+"/A", 1.0 ) ),
      _is_compressible(is_compressible)
  {
    // Force the user to set A
    if( !input.have_variable("Physics/"+physics_name+"/A") )
      {
        std::cerr << "Error: Must specify initial area for "+physics_name << std::endl
                  << "       Input the option Physics/"+physics_name+"/A" << std::endl;
        libmesh_error();
      }

    this->_bc_handler = new SolidMechanicsBCHandling( physics_name, input );

    this->_ic_handler = new GenericICHandler(physics_name, input);

    return;
  }

  template<typename StressStrainLaw>
  ElasticCable<StressStrainLaw>::~ElasticCable()
  {
    return;
  }


  template<typename StressStrainLaw>
  void ElasticCable<StressStrainLaw>::init_variables( libMesh::FEMSystem* system )
  {
    // is_2D = false, is_3D = true
    _disp_vars.init(system,false,true);

    return;
  }


  template<typename StressStrainLaw>
  void ElasticCable<StressStrainLaw>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Tell the system to march temperature forward in time
    system->time_evolving(_disp_vars.u_var());
    system->time_evolving(_disp_vars.v_var());
    system->time_evolving(_disp_vars.w_var());

    return;
  }

  template<typename StressStrainLaw>
  void ElasticCable<StressStrainLaw>::init_context( AssemblyContext& context )
  {
    context.get_element_fe(_disp_vars.u_var())->get_JxW();
    context.get_element_fe(_disp_vars.u_var())->get_phi();
    context.get_element_fe(_disp_vars.u_var())->get_dphidxi();

    // Need for constructing metric tensors
    context.get_element_fe(_disp_vars.u_var())->get_dxyzdxi();
    context.get_element_fe(_disp_vars.u_var())->get_dxidx();
    context.get_element_fe(_disp_vars.v_var())->get_dxidy();
    context.get_element_fe(_disp_vars.w_var())->get_dxidz();

    return;
  }

  template<typename StressStrainLaw>
  void ElasticCable<StressStrainLaw>::register_postprocessing_vars( const GetPot& input,
                                                                    PostProcessedQuantities<libMesh::Real>& postprocessing )
  {
    std::string section = "Physics/"+elastic_cable+"/output_vars";

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
                _stress_indices.resize(3);

                this->_stress_indices[0] = postprocessing.register_quantity("stress_xx");

                this->_stress_indices[1] = postprocessing.register_quantity("stress_xy");

                this->_stress_indices[2] = postprocessing.register_quantity("stress_yy");
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
                std::cerr << "Error: Invalue output_vars value for "+elastic_cable << std::endl
                          << "       Found " << name << std::endl
                          << "       Acceptable values are: stress" << std::endl
                          << "                              strain" << std::endl;
                libmesh_error();
              }
          }
      }
    return;
  }

  template<typename StressStrainLaw>
  void ElasticCable<StressStrainLaw>::element_time_derivative( bool compute_jacobian,
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

    //Grab the Jacobian matrix as submatrices
    //libMesh::DenseMatrix<libMesh::Number> &K = context.get_elem_jacobian();
    libMesh::DenseSubMatrix<libMesh::Number> &Kuu = context.get_elem_jacobian(0,0);
    libMesh::DenseSubMatrix<libMesh::Number> &Kuv = context.get_elem_jacobian(0,1);
    libMesh::DenseSubMatrix<libMesh::Number> &Kuw = context.get_elem_jacobian(0,2);
    libMesh::DenseSubMatrix<libMesh::Number> &Kvu = context.get_elem_jacobian(1,0);
    libMesh::DenseSubMatrix<libMesh::Number> &Kvv = context.get_elem_jacobian(1,1);
    libMesh::DenseSubMatrix<libMesh::Number> &Kvw = context.get_elem_jacobian(1,2);
    libMesh::DenseSubMatrix<libMesh::Number> &Kwu = context.get_elem_jacobian(2,0);
    libMesh::DenseSubMatrix<libMesh::Number> &Kwv = context.get_elem_jacobian(2,1);
    libMesh::DenseSubMatrix<libMesh::Number> &Kww = context.get_elem_jacobian(2,2);


    unsigned int n_qpoints = context.get_element_qrule().n_points();

    // All shape function gradients are w.r.t. master element coordinates
    const std::vector<std::vector<libMesh::Real> >& dphi_dxi = context.get_element_fe(_disp_vars.u_var())->get_dphidxi();

    const libMesh::DenseSubVector<libMesh::Number>& u_coeffs = context.get_elem_solution( _disp_vars.u_var() );
    const libMesh::DenseSubVector<libMesh::Number>& v_coeffs = context.get_elem_solution( _disp_vars.v_var() );
    const libMesh::DenseSubVector<libMesh::Number>& w_coeffs = context.get_elem_solution( _disp_vars.w_var() );

    // Need these to build up the covariant and contravariant metric tensors
    const std::vector<libMesh::RealGradient>& dxdxi  = context.get_element_fe(_disp_vars.u_var())->get_dxyzdxi();

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

        this->compute_metric_tensors( qp, *(context.get_element_fe(_disp_vars.u_var())), context,
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

            for( unsigned int alpha = 0; alpha < dim; alpha++ )
              {
                for( unsigned int beta = 0; beta < dim; beta++ )
                  {
                    Fu(i) -= tau(alpha,beta)*_A*( (grad_x(beta) + grad_u(beta))*u_gradphi(alpha) ) * jac;

                    Fv(i) -= tau(alpha,beta)*_A*( (grad_y(beta) + grad_v(beta))*u_gradphi(alpha) ) * jac;

                    Fw(i) -= tau(alpha,beta)*_A*( (grad_z(beta) + grad_w(beta))*u_gradphi(alpha) ) * jac;
                  }
              }
          }

        if( compute_jacobian )
          {
            for(unsigned int i=0; i != n_u_dofs; i++)
              {
                libMesh::RealGradient u_gradphi_I( dphi_dxi[i][qp] );
                for(unsigned int j=0; j != n_u_dofs; j++)
                  {
                    libMesh::RealGradient u_gradphi_J( dphi_dxi[j][qp] );

                    const libMesh::Real diag_term = _A*jac*tau(0,0)*( u_gradphi_J(0)*u_gradphi_I(0));

                    Kuu(i,j) -= diag_term;

                    Kvv(i,j) -= diag_term;

                    Kww(i,j) -= diag_term;

                    const libMesh::Real dgamma_du = ( u_gradphi_J(0)*(grad_x(0)+grad_u(0)) );

                    const libMesh::Real dgamma_dv = ( u_gradphi_J(0)*(grad_y(0)+grad_v(0)) );

                    const libMesh::Real dgamma_dw = ( u_gradphi_J(0)*(grad_z(0)+grad_w(0)) );

                    const libMesh::Real C1 = _A*jac*C(0,0,0,0);

                    const libMesh::Real x_term = C1*( (grad_x(0)+grad_u(0))*u_gradphi_I(0) );

                    const libMesh::Real y_term = C1*( (grad_y(0)+grad_v(0))*u_gradphi_I(0) );

                    const libMesh::Real z_term = C1*( (grad_z(0)+grad_w(0))*u_gradphi_I(0) );

                    Kuu(i,j) -= x_term*dgamma_du;

                    Kuv(i,j) -= x_term*dgamma_dv;

                    Kuw(i,j) -= x_term*dgamma_dw;

                    Kvu(i,j) -= y_term*dgamma_du;

                    Kvv(i,j) -= y_term*dgamma_dv;

                    Kvw(i,j) -= y_term*dgamma_dw;

                    Kwu(i,j) -= z_term*dgamma_du;

                    Kwv(i,j) -= z_term*dgamma_dv;

                    Kww(i,j) -= z_term*dgamma_dw;

                  }
              }
            //libmesh_not_implemented();
          }
      }

    return;
  }

  template<typename StressStrainLaw>
  void ElasticCable<StressStrainLaw>::side_time_derivative( bool compute_jacobian,
                                                            AssemblyContext& context,
                                                            CachedValues& cache )
  {
    std::vector<BoundaryID> ids = context.side_boundary_ids();

    for( std::vector<BoundaryID>::const_iterator it = ids.begin();
         it != ids.end(); it++ )
      {
        libmesh_assert (*it != libMesh::BoundaryInfo::invalid_id);

        _bc_handler->apply_neumann_bcs( context, cache, compute_jacobian, *it );
      }

    return;
  }

  template<typename StressStrainLaw>
  void ElasticCable<StressStrainLaw>::mass_residual( bool /*compute_jacobian*/,
                                                     AssemblyContext& /*context*/,
                                                     CachedValues& /*cache*/ )
  {
    libmesh_not_implemented();
    return;
  }

  template<typename StressStrainLaw>
  void ElasticCable<StressStrainLaw>::compute_postprocessed_quantity( unsigned int quantity_index,
                                                                      const AssemblyContext& context,
                                                                      const libMesh::Point& point,
                                                                      libMesh::Real& value )
  {
    bool is_strain = false;
    if( !_strain_indices.empty() )
      is_strain = ( _strain_indices[0] == quantity_index );

    if( is_strain )
      {
        const unsigned int n_u_dofs = context.get_dof_indices(_disp_vars.u_var()).size();

        const libMesh::DenseSubVector<libMesh::Number>& u_coeffs = context.get_elem_solution( _disp_vars.u_var() );
        const libMesh::DenseSubVector<libMesh::Number>& v_coeffs = context.get_elem_solution( _disp_vars.v_var() );
        const libMesh::DenseSubVector<libMesh::Number>& w_coeffs = context.get_elem_solution( _disp_vars.w_var() );

        // Build new FE for the current point. We need this to build tensors at point.
        libMesh::AutoPtr<libMesh::FEGenericBase<libMesh::Real> > fe_new =  this->build_new_fe( context.get_elem(), context.get_element_fe(_disp_vars.u_var()), point );

        const std::vector<std::vector<libMesh::Real> >& dphi_dxi =  fe_new->get_dphidxi();

        // Need these to build up the covariant and contravariant metric tensors
        const std::vector<libMesh::RealGradient>& dxdxi  = fe_new->get_dxyzdxi();

        // Gradients are w.r.t. master element coordinates
        libMesh::Gradient grad_u, grad_v, grad_w;
        for( unsigned int d = 0; d < n_u_dofs; d++ )
          {
            libMesh::RealGradient u_gradphi( dphi_dxi[d][0] );
            grad_u += u_coeffs(d)*u_gradphi;
            grad_v += v_coeffs(d)*u_gradphi;
            grad_w += w_coeffs(d)*u_gradphi;
          }

        libMesh::RealGradient grad_x( dxdxi[0](0) );
        libMesh::RealGradient grad_y( dxdxi[0](1) );
        libMesh::RealGradient grad_z( dxdxi[0](2) );

        libMesh::TensorValue<libMesh::Real> a_cov, a_contra, A_cov, A_contra;
        libMesh::Real lambda_sq = 0;

        this->compute_metric_tensors(0, *fe_new, context, grad_u, grad_v, grad_w, a_cov, a_contra, A_cov, A_contra, lambda_sq );

        // We have everything we need for strain now, so check if we are computing strain
        if( is_strain )
          {
            if( _strain_indices[0] == quantity_index )
              {
                value = 0.5*(A_cov(0,0) - a_cov(0,0));
              }
            else
              {
                //Wat?!
                libmesh_error();
              }
            return;
          }
      }

    return;
  }

  template<typename StressStrainLaw>
  libMesh::AutoPtr<libMesh::FEGenericBase<libMesh::Real> > ElasticCable<StressStrainLaw>::build_new_fe( const libMesh::Elem& elem,
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

  template<typename StressStrainLaw>
  void ElasticCable<StressStrainLaw>::compute_metric_tensors( unsigned int qp,
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
    a_contra(0,0) = dxi*dxi;
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

    return;
  }

} // end namespace GRINS
