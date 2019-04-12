//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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
#include "grins/vorticity.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/assembly_context.h"
#include "grins/variables_parsing.h"
#include "grins/variable_warehouse.h"
#include "grins/multi_component_vector_variable.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"
#include "libmesh/fe_base.h"
#include "libmesh/elem.h"

namespace GRINS
{
  Vorticity::Vorticity( const std::string& qoi_name )
    : QoIBase(qoi_name)
  {
    return;
  }

  Vorticity::~Vorticity()
  {
    return;
  }

  QoIBase* Vorticity::clone() const
  {
    return new Vorticity( *this );
  }

  void Vorticity::init
  (const GetPot& input,
   const MultiphysicsSystem& /*system*/,
   unsigned int /*qoi_num*/ )
  {
    // Extract subdomain on which to compute to qoi
    int num_ids = input.vector_variable_size( "QoI/Vorticity/enabled_subdomains" );

    if( num_ids == 0 )
      {
        std::cerr << "Error: Must specify at least one subdomain id on which to compute vorticity." << std::endl;
        libmesh_error();
      }

    for( int i = 0; i < num_ids; i++ )
      {
        libMesh::subdomain_id_type s_id = input( "QoI/Vorticity/enabled_subdomains", -1, i );
        _subdomain_ids.insert( s_id );
      }

    _flow_vars = &(GRINSPrivate::VariableWarehouse::get_variable_subclass<VelocityVariable>(VariablesParsing::velocity_variable_name(input,"Vorticity",VariablesParsing::QOI)));
  }

  void Vorticity::init_context( AssemblyContext& context )
  {
    libMesh::FEBase* u_fe = NULL;
    libMesh::FEBase* v_fe = NULL;

    context.get_element_fe<libMesh::Real>(_flow_vars->u(), u_fe);
    context.get_element_fe<libMesh::Real>(_flow_vars->v(), v_fe);

    u_fe->get_dphi();
    u_fe->get_JxW();

    v_fe->get_dphi();

    return;
  }

  void Vorticity::element_qoi( AssemblyContext& context,
                               const unsigned int qoi_index )
  {

    if( _subdomain_ids.find( (&context.get_elem())->subdomain_id() ) != _subdomain_ids.end() )
      {
        libMesh::FEBase* element_fe;
        context.get_element_fe<libMesh::Real>(_flow_vars->u(), element_fe);
        const std::vector<libMesh::Real> &JxW = element_fe->get_JxW();

        unsigned int n_qpoints = context.get_element_qrule().n_points();

        /*! \todo Need to generalize this to the multiple QoI case */
        libMesh::Number& qoi = context.get_qois()[qoi_index];

        for( unsigned int qp = 0; qp != n_qpoints; qp++ )
          {
            libMesh::Gradient grad_u = 0.;
            libMesh::Gradient grad_v = 0.;
            context.interior_gradient( _flow_vars->u(), qp, grad_u );
            context.interior_gradient( _flow_vars->v(), qp, grad_v );
            qoi += (grad_v(0) - grad_u(1)) * JxW[qp];
          }
      }

    return;
  }

  void Vorticity::element_qoi_derivative( AssemblyContext& context,
                                          const unsigned int qoi_index )
  {

    if( _subdomain_ids.find( (&context.get_elem())->subdomain_id() ) != _subdomain_ids.end() )
      {
        // Element
        libMesh::FEBase* element_fe;
        context.get_element_fe<libMesh::Real>(_flow_vars->u(), element_fe);

        // Jacobian times integration weights
        const std::vector<libMesh::Real> &JxW = element_fe->get_JxW();

        // Grad of basis functions
        const std::vector<std::vector<libMesh::RealGradient> >& du_phi =
          context.get_element_fe(_flow_vars->u())->get_dphi();
        const std::vector<std::vector<libMesh::RealGradient> >& dv_phi =
          context.get_element_fe(_flow_vars->v())->get_dphi();

        // Local DOF count and quadrature point count
        const unsigned int n_T_dofs = context.get_dof_indices(0).size();
        unsigned int n_qpoints = context.get_element_qrule().n_points();

        // Warning: we assume here that vorticity is the only QoI!
        // This should be consistent with the assertion in grins_mesh_adaptive_solver.C
        /*! \todo Need to generalize this to the multiple QoI case */
        libMesh::DenseSubVector<libMesh::Number> &Qu = context.get_qoi_derivatives(qoi_index, _flow_vars->u());
        libMesh::DenseSubVector<libMesh::Number> &Qv = context.get_qoi_derivatives(qoi_index, _flow_vars->v());

        // Integration loop
        for( unsigned int qp = 0; qp != n_qpoints; qp++ )
          {
            for( unsigned int i = 0; i != n_T_dofs; i++ )
              {
                Qu(i) += - dv_phi[i][qp](1) * JxW[qp];
                Qv(i) += du_phi[i][qp](0) * JxW[qp];
              }
          }
      }

    return;
  }

} //namespace GRINS
