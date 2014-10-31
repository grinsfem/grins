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

// libMesh
// libMesh
#include "libmesh/getpot.h"
#include "libmesh/quadrature.h"
#include "libmesh/boundary_info.h"

namespace GRINS
{
  template<typename ElasticityTensor>
  ElasticMembrane<ElasticityTensor>::ElasticMembrane( const GRINS::PhysicsName& physics_name, const GetPot& input )
    : Physics(physics_name,input)
  {
    return;
  }
  
  template<typename ElasticityTensor>
  ElasticMembrane<ElasticityTensor>::~ElasticMembrane()
  {
    return;
  }

  template<typename ElasticityTensor>
  void ElasticMembrane<ElasticityTensor>::init_variables( libMesh::FEMSystem* system )
  {
    // This will be the manifold dimension (2), not the spatial dimension (3).
    this->_dim = system->get_mesh().mesh_dimension();

    // is_2D = false, is_3D = true
    _disp_vars.init(system,false,true);

    return;
  }

  template<typename ElasticityTensor>
  void ElasticMembrane<ElasticityTensor>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Tell the system to march temperature forward in time
    system->time_evolving(_disp_vars.u_var());
    system->time_evolving(_disp_vars.v_var());
    system->time_evolving(_disp_vars.w_var());

    return;
  }

  void HeatTransferBase::init_context( AssemblyContext& context )
  {
    context.get_element_fe(_disp_vars.u_var())->get_JxW();
    context.get_element_fe(_disp_vars.u_var())->get_dphi();

    // Need for constructing metric tensors
    context.get_element_fe(_disp_vars.u_var())->get_dxyzdxi();
    context.get_element_fe(_disp_vars.u_var())->get_dxyzdeta();
    context.get_element_fe(_disp_vars.u_var())->get_dxidx();
    context.get_element_fe(_disp_vars.u_var())->get_dxidy();
    context.get_element_fe(_disp_vars.u_var())->get_dxidz();
    context.get_element_fe(_disp_vars.u_var())->get_detadx();
    context.get_element_fe(_disp_vars.u_var())->get_detady();
    context.get_element_fe(_disp_vars.u_var())->get_detadz();
    

    return;
  }

  template<typename ElasticityTensor>
  void ElasticMembrane<ElasticityTensor>::element_time_derivative( bool compute_jacobian,
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

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Gradient grad_u = context.interior_gradient(_disp_vars.u_var(), qp);
        libMesh::Gradient grad_v = context.interior_gradient(_disp_vars.v_var(), qp);
        libMesh::Gradient grad_w = context.interior_gradient(_disp_vars.v_var(), qp);

        libMesh::RealGradient grad_x( dxdxi[qp](0), dxdeta[qp](0) );
        libMesh::RealGradient grad_y( dxdxi[qp](1), dxdeta[qp](1) );
        libMesh::RealGradient grad_z( dxdxi[qp](2), dxdeta[qp](2) );

        libMesh::RealGradient dudxi( grad_u(0), grad_v(0), grad_w(0) );
        libMesh::RealGradient dudeta( grad_u(1), grad_v(1), grad_w(1) );
        
        libMesh::RealGradient dxi( dxidx[qp], dxidy[qp], dxidz[qp] );
        libMesh::RealGradient deta( detadx[qp], detady[qp], detadz[qp] );

        // Covariant metric tensor of reference configuration
        libMesh::TensorValue<libMesh::Real> a_cov( dxdxi[qp]*dxdxi[qp], dxdxi[qp]*dxdeta[qp], 0.0,
                                                   dxeta[qp]*dxdxi[qp], dxdeta[qp]*dxdeta[qp] );

        // Covariant metric tensor of current configuration
        libMesh::TensorValue<libMesh::Real> A_cov( (dxdxi[qp] + dudxi)*(dxdxi[qp] + dudxi),
                                                   (dxdxi[qp] + dudxi)*(dxdeta[qp] + dudeta), 0.0,
                                                   (dxeta[qp] + dudeta)*(dxdxi[qp] + dudxi),
                                                   (dxdeta[qp] + dudeta)*(dxdeta[qp] + dudeta) );

        
        // Contravariant metric tensor of reference configuration
        libMesh::TensorValue<libMesh::Real> a_contra( dxi*dxi, dxi*deta, 0.0
                                                      deta*dxi, deta*deta ); 

        // Strain tensor
        libMesh::TensorValue<libMesh::Real> strain = 0.5*(A_cov - a_cov);

        // Compute stress tensor
        libMesh::TensorValue<libMesh::Real> tau = this->compute_stress(a_contra,strain);

        libMesh::Real jac = JxW[qp];

        for (unsigned int i=0; i != n_u_dofs; i++)
	  {
            for( unsigned int alpha = 0; alpha < 2; alpha++ )
              {
                for( unsigned int beta = 0; beta < 2; beta++ )
                  {
                    Fu(i) -= 0.5*tau(alpha,beta)*( (grad_x(beta) + grad_u(beta))*u_gradphi[i][qp](alpha) +
                                                   (grad_x(alpha) + grad_u(alpha))*u_gradphi[i][qp](beta) )*jac;

                    Fv(i) -= 0.5*tau(alpha,beta)*( (grad_y(beta) + grad_v(beta))*u_gradphi[i][qp](alpha) +
                                                   (grad_y(alpha) + grad_v(alpha))*u_gradphi[i][qp](beta) )*jac;

                    Fw(i) -= 0.5*tau(alpha,beta)*( (grad_z(beta) + grad_w(beta))*u_gradphi[i][qp](alpha) +
                                                   (grad_z(alpha) + grad_w(alpha))*u_gradphi[i][qp](beta) )*jac;
                  }
              }
          }

      }

    return;
  }

  template<typename ElasticityTensor>
  void ElasticMembrane<ElasticityTensor>::side_time_derivative( bool compute_jacobian,
                                              AssemblyContext& context,
                                              CachedValues& cache )
  {
    libmesh_not_implemented();
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

  template<typename ElasticityTensor>
  void ElasticMembrane<ElasticityTensor>::mass_residual( bool compute_jacobian,
                                       AssemblyContext& context,
                                       CachedValues& cache )
  {
    libmesh_not_implemented();
    return;
  }

  template<typename ElasticityTensor>
  libMesh::TensorValue<libMesh::Real> ElasticMembrane<ElasticityTensor>::compute_stress(libMesh::TensorValue<libMesh::Real>& a_contra,
                                                                                        libMesh::TensorValue<libMesh::Real>& strain)
  {
    libMesh::TensorValue<libMesh::Real> tau;

    (this->C).recompute_elasticity(a_contra);

    for( unsigned int alpha = 0; alpha < 2; alpha++ )
      {
        for( unsigned int beta = 0; beta < 2; beta++ )
          {
            for( unsigned int gamma = 0; gamma < 2; gamma++ )
              {
                for( unsigned int delta = 0; delta < 2; delta++ )
                  {
                    tau(alpha,beta) += (this->C(alpha,beta,gamma,delta))*strain(gamma,delta)
                  }
              }
          }
      }

    return tau;
  }

} // end namespace GRINS
