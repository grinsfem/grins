//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "grins/vorticity.h"

// GRINS
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_context.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace GRINS
{
  Vorticity::Vorticity( const GetPot& input )
    : QoIBase()
  {
    this->assemble_qoi_elements = true;
    this->read_input_options( input );

    return;
  }

  Vorticity::~Vorticity()
  {
    return;
  }

  libMesh::AutoPtr<libMesh::DifferentiableQoI> Vorticity::clone()
  {
    return libMesh::AutoPtr<libMesh::DifferentiableQoI>( new Vorticity( *this ) );
  }

  void Vorticity::init( const GetPot& input, const MultiphysicsSystem& system )
  {
    // Grab velocity variable indices
    std::string u_var_name = input("Physics/VariableNames/u_velocity", u_var_name_default);
    std::string v_var_name = input("Physics/VariableNames/v_velocity", v_var_name_default);
    this->_u_var = system.variable_number(u_var_name);
    this->_v_var = system.variable_number(v_var_name);

    return;
  }

  void Vorticity::read_input_options( const GetPot& input )
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

    return;
  }

  void Vorticity::element_qoi( libMesh::DiffContext& context, const libMesh::QoISet& )
  {
    libMesh::FEMContext &c = libmesh_cast_ref<libMesh::FEMContext&>(context);

    if( _subdomain_ids.find( (c.elem)->subdomain_id() ) != _subdomain_ids.end() )
      {
	libMesh::FEBase* element_fe;
	c.get_element_fe<libMesh::Real>(this->_u_var, element_fe);
	const std::vector<libMesh::Real> &JxW = element_fe->get_JxW();

	unsigned int n_qpoints = (c.get_element_qrule())->n_points();

	/*! \todo Need to generalize this to the multiple QoI case */
	libMesh::Number& qoi = c.elem_qoi[0];

	for( unsigned int qp = 0; qp != n_qpoints; qp++ )
	  {
	    libMesh::Gradient grad_u = 0.;
	    libMesh::Gradient grad_v = 0.;
	    c.interior_gradient( this->_u_var, qp, grad_u );
	    c.interior_gradient( this->_v_var, qp, grad_v );
	    qoi += (grad_v(0) - grad_u(1)) * JxW[qp];
	  }
      }

    return;
  }

  void Vorticity::element_qoi_derivative( libMesh::DiffContext &context, const libMesh::QoISet & )
  {
    libMesh::FEMContext &c = libmesh_cast_ref<libMesh::FEMContext&>(context);

    if( _subdomain_ids.find( (c.elem)->subdomain_id() ) != _subdomain_ids.end() )
      {
	// Element
	libMesh::FEBase* element_fe;
	c.get_element_fe<libMesh::Real>(this->_u_var, element_fe);

	// Jacobian times integration weights
	const std::vector<libMesh::Real> &JxW = element_fe->get_JxW();

	// Grad of basis functions
	const std::vector<std::vector<libMesh::RealGradient> >& du_phi =
	  c.element_fe_var[_u_var]->get_dphi();
	const std::vector<std::vector<libMesh::RealGradient> >& dv_phi =
	  c.element_fe_var[_v_var]->get_dphi();

	// Local DOF count and quadrature point count
	const unsigned int n_T_dofs = c.dof_indices_var[0].size();
	unsigned int n_qpoints = (c.get_element_qrule())->n_points();  

	// Warning: we assume here that vorticity is the only QoI!
	// This should be consistent with the assertion in grins_mesh_adaptive_solver.C
	/*! \todo Need to generalize this to the multiple QoI case */
	libMesh::DenseSubVector<Number> &Qu = *c.elem_qoi_subderivatives[0][0];
	libMesh::DenseSubVector<Number> &Qv = *c.elem_qoi_subderivatives[0][1];

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
