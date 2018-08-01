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
#include "grins/elastic_cable.h"

// GRINS
#include "grins_config.h"
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
  ElasticCable<StressStrainLaw>::ElasticCable( const PhysicsName& physics_name, const GetPot& input,
                                               bool is_compressible )
    : ElasticCableBase<StressStrainLaw>(physics_name,input,is_compressible)
  {
    this->_ic_handler = new GenericICHandler(physics_name, input);
  }

  template<typename StressStrainLaw>
  void ElasticCable<StressStrainLaw>::register_postprocessing_vars( const GetPot& input,
                                                                    PostProcessedQuantities<libMesh::Real>& postprocessing )
  {
    std::string section = "Physics/"+PhysicsNaming::elastic_cable()+"/output_vars";

    if( input.have_variable(section) )
      {
        unsigned int n_vars = input.vector_variable_size(section);

        for( unsigned int v = 0; v < n_vars; v++ )
          {
            std::string name = input(section,"DIE!",v);

            if( name == std::string("stress") )
              {
                // sigma_xx only, i.e., locally normal to the cutting plane
                // sigma_yy=sigma_zz = 0 by assumption of this Physics
                _stress_indices.resize(1);

                this->_stress_indices[0] = postprocessing.register_quantity("cable_stress");

              }
            else if( name == std::string("strain") )
              {
                // eps_xx only, i.e., locally normal to the cutting plane
                // eps_yy=eps_zz = 0 by assumption of this Physics
                _strain_indices.resize(1);

                this->_strain_indices[0] = postprocessing.register_quantity("cable_strain");

              }
            else if( name == std::string("force") )
              {
                // force_x only, i.e., locally normal to the cutting plane
                // force_y=force_z=0 by assumption of this Physics
                _force_indices.resize(1);

                this->_force_indices[0] = postprocessing.register_quantity("cable_force");

              }
            else
              {
                std::cerr << "Error: Invalue output_vars value for "+PhysicsNaming::elastic_cable() << std::endl
                          << "       Found " << name << std::endl
                          << "       Acceptable values are: stress" << std::endl
                          << "                              strain" << std::endl
                          << "                              force " << std::endl;
                libmesh_error();
              }
          }
      }
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

    bool is_stress = false;
    if( !_stress_indices.empty() )
      is_stress = ( _stress_indices[0] == quantity_index );

    bool is_force = false;
    if( !_force_indices.empty() )
      is_force = ( _force_indices[0] == quantity_index );

    if( is_strain || is_stress || is_force )
      {
        const unsigned int n_u_dofs = context.get_dof_indices(this->_disp_vars.u()).size();

        const libMesh::DenseSubVector<libMesh::Number>& u_coeffs = context.get_elem_solution( this->_disp_vars.u() );
        const libMesh::DenseSubVector<libMesh::Number>* v_coeffs = NULL;

        const libMesh::DenseSubVector<libMesh::Number>* w_coeffs = NULL;

        if( this->_disp_vars.dim() >= 2 )
          v_coeffs = &context.get_elem_solution( this->_disp_vars.v() );

        if( this->_disp_vars.dim() == 3 )
          w_coeffs = &context.get_elem_solution( this->_disp_vars.w() );

        // Build new FE for the current point. We need this to build tensors at point.
        std::unique_ptr<libMesh::FEGenericBase<libMesh::Real> > fe_new =
          this->build_new_fe( &context.get_elem(), this->get_fe(context), point );

        const std::vector<std::vector<libMesh::Real> >& dphi_dxi =  fe_new->get_dphidxi();

        // Need these to build up the covariant and contravariant metric tensors
        const std::vector<libMesh::RealGradient>& dxdxi  = fe_new->get_dxyzdxi();

        // Gradients are w.r.t. master element coordinates
        libMesh::Gradient grad_u, grad_v, grad_w;
        for( unsigned int d = 0; d < n_u_dofs; d++ )
          {
            libMesh::RealGradient u_gradphi( dphi_dxi[d][0] );
            grad_u += u_coeffs(d)*u_gradphi;

            if( this->_disp_vars.dim() >= 2 )
              grad_v += (*v_coeffs)(d)*u_gradphi;

            if( this->_disp_vars.dim() == 3 )
              grad_w += (*w_coeffs)(d)*u_gradphi;
          }

        libMesh::RealGradient grad_x( dxdxi[0](0) );
        libMesh::RealGradient grad_y( dxdxi[0](1) );
        libMesh::RealGradient grad_z( dxdxi[0](2) );

        libMesh::TensorValue<libMesh::Real> a_cov, a_contra, A_cov, A_contra;
        libMesh::Real lambda_sq = 0;

        this->compute_metric_tensors(0, *fe_new, context,
                                     grad_u, grad_v, grad_w,
                                     a_cov, a_contra, A_cov, A_contra, lambda_sq );

        libMesh::Real det_a = a_cov(0,0)*a_cov(1,1) - a_cov(0,1)*a_cov(1,0);
        libMesh::Real det_A = A_cov(0,0)*A_cov(1,1) - A_cov(0,1)*A_cov(1,0);

        libMesh::Real I3 = lambda_sq*det_A/det_a;

        libMesh::TensorValue<libMesh::Real> tau;
        this->_stress_strain_law.compute_stress(1,a_contra,a_cov,A_contra,A_cov,tau);

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

        if( is_stress )
          {
            if( _stress_indices[0] == quantity_index )
              {
                // Need to convert to Cauchy stress
                value = tau(0,0)/std::sqrt(I3);
              }
            else
              {
                libmesh_error();
              }
          }
        if( is_force )
          {

            if( _force_indices[0] == quantity_index )
              {
                //This force is in deformed configuration and will be significantly influenced by large strains
                value = tau(0,0)/std::sqrt(I3)*this->_A;
              }
            else
              {
                libmesh_error();
              }
          }
      }
  }

  template<typename StressStrainLaw>
  void ElasticCable<StressStrainLaw>::element_time_derivative
  ( bool compute_jacobian, AssemblyContext & context )
  {
    unsigned int u_var = this->_disp_vars.u();

    const unsigned int n_u_dofs = context.get_dof_indices(u_var).size();

    const std::vector<libMesh::Real> &JxW =
      this->get_fe(context)->get_JxW();

    const MultiphysicsSystem & system = context.get_multiphysics_system();

    unsigned int u_dot_var = system.get_second_order_dot_var(u_var);

    // Residuals that we're populating
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(u_dot_var);
    libMesh::DenseSubVector<libMesh::Number>* Fv = NULL;
    libMesh::DenseSubVector<libMesh::Number>* Fw = NULL;

    libMesh::DenseSubMatrix<libMesh::Number>& Kuu = context.get_elem_jacobian(u_dot_var,u_var);
    libMesh::DenseSubMatrix<libMesh::Number>* Kuv = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kuw = NULL;

    libMesh::DenseSubMatrix<libMesh::Number>* Kvu = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kvv = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kvw = NULL;

    libMesh::DenseSubMatrix<libMesh::Number>* Kwu = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kwv = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kww = NULL;

    unsigned int v_var = libMesh::invalid_uint;
    unsigned int v_dot_var = libMesh::invalid_uint;
    if( this->_disp_vars.dim() >= 2 )
      {
        v_var = this->_disp_vars.v();
        v_dot_var = system.get_second_order_dot_var(v_var);

        Fv = &context.get_elem_residual(v_dot_var);

        Kuv = &context.get_elem_jacobian(u_dot_var,v_var);
        Kvu = &context.get_elem_jacobian(v_dot_var,u_var);
        Kvv = &context.get_elem_jacobian(v_dot_var,v_var);
      }

    unsigned int w_var = libMesh::invalid_uint;
    unsigned int w_dot_var = libMesh::invalid_uint;
    if( this->_disp_vars.dim() == 3 )
      {
        w_var = this->_disp_vars.w();
        w_dot_var = system.get_second_order_dot_var(w_var);

        Fw = &context.get_elem_residual(this->_disp_vars.w());

        Kuw = &context.get_elem_jacobian(u_dot_var,w_var);
        Kvw = &context.get_elem_jacobian(v_dot_var,w_var);
        Kwu = &context.get_elem_jacobian(w_dot_var,u_var);
        Kwv = &context.get_elem_jacobian(w_dot_var,v_var);
        Kww = &context.get_elem_jacobian(w_dot_var,w_var);
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    // All shape function gradients are w.r.t. master element coordinates
    const std::vector<std::vector<libMesh::Real> >& dphi_dxi = this->get_fe(context)->get_dphidxi();

    // Need these to build up the covariant and contravariant metric tensors
    const std::vector<libMesh::RealGradient>& dxdxi  = this->get_fe(context)->get_dxyzdxi();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Real jac = JxW[qp];
        libMesh::Gradient grad_u,grad_v,grad_w;
        this->get_grad_disp(context, qp, grad_u,grad_v,grad_w);

        libMesh::TensorValue<libMesh::Real> tau;
        ElasticityTensor C;
        this->get_stress_and_elasticity(context,qp,grad_u,grad_v,grad_w,tau,C);

        libMesh::RealGradient grad_x( dxdxi[qp](0) );
        libMesh::RealGradient grad_y( dxdxi[qp](1) );
        libMesh::RealGradient grad_z( dxdxi[qp](2) );

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            libMesh::RealGradient u_gradphi( dphi_dxi[i][qp] );

            const libMesh::Real res_term = tau(0,0)*this->_A*jac*u_gradphi(0);

            Fu(i) += res_term*(grad_x(0) + grad_u(0));

            if( this->_disp_vars.dim() >= 2 )
              (*Fv)(i) += res_term*(grad_y(0) + grad_v(0));

            if( this->_disp_vars.dim() == 3 )
              (*Fw)(i) += res_term*(grad_z(0) + grad_w(0));
          }

        if( compute_jacobian )
          {
            for(unsigned int i=0; i != n_u_dofs; i++)
              {
                libMesh::RealGradient u_gradphi_I( dphi_dxi[i][qp] );
                for(unsigned int j=0; j != n_u_dofs; j++)
                  {
                    libMesh::RealGradient u_gradphi_J( dphi_dxi[j][qp] );

                    const libMesh::Real diag_term =
                      this->_A*jac*tau(0,0)*( u_gradphi_J(0)*u_gradphi_I(0))*context.get_elem_solution_derivative();

                    Kuu(i,j) += diag_term;

                    if( this->_disp_vars.dim() >= 2 )
                      (*Kvv)(i,j) += diag_term;

                    if( this->_disp_vars.dim() == 3 )
                      (*Kww)(i,j) += diag_term;

                    const libMesh::Real dgamma_du = ( u_gradphi_J(0)*(grad_x(0)+grad_u(0)) );

                    const libMesh::Real C1 = this->_A*jac*C(0,0,0,0)*context.get_elem_solution_derivative();

                    const libMesh::Real x_term = C1*( (grad_x(0)+grad_u(0))*u_gradphi_I(0) );

                    Kuu(i,j) += x_term*dgamma_du;

                    libMesh::Real y_term = 0.0;
                    libMesh::Real dgamma_dv = 0.0;

                    if( this->_disp_vars.dim() >= 2 )
                      {
                        dgamma_dv = ( u_gradphi_J(0)*(grad_y(0)+grad_v(0)) );
                        y_term = C1*( (grad_y(0)+grad_v(0))*u_gradphi_I(0) );

                        (*Kuv)(i,j) += x_term*dgamma_dv;
                        (*Kvu)(i,j) += y_term*dgamma_du;
                        (*Kvv)(i,j) += y_term*dgamma_dv;
                      }

                    libMesh::Real z_term = 0.0;

                    if( this->_disp_vars.dim() == 3 )
                      {
                        const libMesh::Real dgamma_dw = ( u_gradphi_J(0)*(grad_z(0)+grad_w(0)) );
                        z_term = C1*( (grad_z(0)+grad_w(0))*u_gradphi_I(0) );

                        (*Kuw)(i,j) += x_term*dgamma_dw;
                        (*Kvw)(i,j) += y_term*dgamma_dw;
                        (*Kwu)(i,j) += z_term*dgamma_du;
                        (*Kwv)(i,j) += z_term*dgamma_dv;
                        (*Kww)(i,j) += z_term*dgamma_dw;
                      }

                  } // end j-loop
              } // end i-loop
          } // end if(compute_jacobian)
      } // end qp loop
  }

  template<typename StressStrainLaw>
  void ElasticCable<StressStrainLaw>::get_grad_disp( const AssemblyContext & context,
                                                     unsigned int qp,
                                                     libMesh::Gradient & grad_u,
                                                     libMesh::Gradient & grad_v,
                                                     libMesh::Gradient & grad_w )
  {
    const int n_u_dofs = context.get_dof_indices(this->_disp_vars.u()).size();

    // All shape function gradients are w.r.t. master element coordinates
    const std::vector<std::vector<libMesh::Real> >& dphi_dxi = this->get_fe(context)->get_dphidxi();

    const libMesh::DenseSubVector<libMesh::Number>& u_coeffs = context.get_elem_solution( this->_disp_vars.u() );
    const libMesh::DenseSubVector<libMesh::Number>* v_coeffs = NULL;
    const libMesh::DenseSubVector<libMesh::Number>* w_coeffs = NULL;

    if( this->_disp_vars.dim() >= 2 )
      v_coeffs = &context.get_elem_solution( this->_disp_vars.v() );

    if( this->_disp_vars.dim() == 3 )
      w_coeffs = &context.get_elem_solution( this->_disp_vars.w() );

    // Compute gradients  w.r.t. master element coordinates
    for( int d = 0; d < n_u_dofs; d++ )
      {
        libMesh::RealGradient u_gradphi( dphi_dxi[d][qp] );
        grad_u += u_coeffs(d)*u_gradphi;

        if( this->_disp_vars.dim() >= 2 )
          grad_v += (*v_coeffs)(d)*u_gradphi;

        if( this->_disp_vars.dim() == 3 )
          grad_w += (*w_coeffs)(d)*u_gradphi;
      }
  }


  template<typename StressStrainLaw>
  void ElasticCable<StressStrainLaw>::get_stress_and_elasticity( const AssemblyContext & context,
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
    const unsigned int dim = 1; // The cable dimension is always 1 for this physics
    this->_stress_strain_law.compute_stress_and_elasticity(dim,a_contra,a_cov,A_contra,A_cov,tau,C);
  }

} // end namespace GRINS
