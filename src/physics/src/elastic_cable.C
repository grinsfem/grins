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
#include "grins/solid_mechanics_bc_handling.h"
#include "grins/generic_ic_handler.h"
#include "grins/postprocessed_quantities.h"

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
    this->_bc_handler = new SolidMechanicsBCHandling( physics_name, input );

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
  void ElasticCable<StressStrainLaw>::side_time_derivative( bool compute_jacobian,
                                                            AssemblyContext& context,
                                                            CachedValues& cache )
  {
    std::vector<BoundaryID> ids = context.side_boundary_ids();

    for( std::vector<BoundaryID>::const_iterator it = ids.begin();
	 it != ids.end(); it++ )
      {
	libmesh_assert (*it != libMesh::BoundaryInfo::invalid_id);

	this->_bc_handler->apply_neumann_bcs( context, cache, compute_jacobian, *it );
      }
  }


  template<typename StressStrainLaw>
  void ElasticCable<StressStrainLaw>::mass_residual( bool compute_jacobian,
                                                     AssemblyContext& context,
                                                     CachedValues& /*cache*/ )
  {
    const unsigned int n_u_dofs = context.get_dof_indices(this->_disp_vars.u()).size();

    const std::vector<libMesh::Real> &JxW =
      this->get_fe(context)->get_JxW();

    const std::vector<std::vector<libMesh::Real> >& u_phi =
      this->get_fe(context)->get_phi();

    // Residuals that we're populating
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_disp_vars.u());
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_disp_vars.v());
    libMesh::DenseSubVector<libMesh::Number> &Fw = context.get_elem_residual(this->_disp_vars.w());

    libMesh::DenseSubMatrix<libMesh::Number>& Kuu = context.get_elem_jacobian(this->_disp_vars.u(),this->_disp_vars.u());
    libMesh::DenseSubMatrix<libMesh::Number>& Kvv = context.get_elem_jacobian(this->_disp_vars.v(),this->_disp_vars.v());
    libMesh::DenseSubMatrix<libMesh::Number>& Kww = context.get_elem_jacobian(this->_disp_vars.w(),this->_disp_vars.w());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Real jac = JxW[qp];

        libMesh::Real u_ddot, v_ddot, w_ddot;
        context.interior_accel( this->_disp_vars.u(), qp, u_ddot );
        context.interior_accel( this->_disp_vars.v(), qp, v_ddot );
        context.interior_accel( this->_disp_vars.w(), qp, w_ddot );

        for (unsigned int i=0; i != n_u_dofs; i++)
	  {
            Fu(i) += this->_rho*this->_A*u_ddot*u_phi[i][qp]*jac;
            Fv(i) += this->_rho*this->_A*v_ddot*u_phi[i][qp]*jac;
            Fw(i) += this->_rho*this->_A*w_ddot*u_phi[i][qp]*jac;

            if( compute_jacobian )
              {
                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    libMesh::Real jac_term = this->_rho*this->_A*u_phi[i][qp]*u_phi[j][qp]*jac;
                    jac_term *= context.get_elem_solution_accel_derivative();

                    Kuu(i,j) += jac_term;
                    Kvv(i,j) += jac_term;
                    Kww(i,j) += jac_term;
                  }
              }
          }
      }
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
        const libMesh::DenseSubVector<libMesh::Number>& v_coeffs = context.get_elem_solution( this->_disp_vars.v() );
        const libMesh::DenseSubVector<libMesh::Number>& w_coeffs = context.get_elem_solution( this->_disp_vars.w() );

        // Build new FE for the current point. We need this to build tensors at point.
        libMesh::AutoPtr<libMesh::FEGenericBase<libMesh::Real> > fe_new =  this->build_new_fe( context.get_elem(), this->get_fe(context), point );

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

} // end namespace GRINS
