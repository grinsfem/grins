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
#include "grins/elastic_cable_constant_gravity.h"

// GRINS
#include "grins_config.h"
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/quadrature.h"
#include "libmesh/boundary_info.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  ElasticCableConstantGravity::ElasticCableConstantGravity( const GRINS::PhysicsName& physics_name, const GetPot& input )
  	  : Physics(physics_name,input),
  	    _disp_vars(input,physics_name),
		_A( input("Physics/"+physics_name+"/A", 1.0 ) ),
		_rho( input("Physics/"+physics_name+"/rho", 1.0 ) )
  {
	if( !input.have_variable("Physics/ElasticCableConstantGravity/A") )
	{
		std::cerr << "Error: Must input area for ElasticCableConstantGravity." << std::endl
				  << "       Please set Physics/ElasticCableConstantGravity/A." << std::endl;
		libmesh_error();
	}
    if( !input.have_variable("Physics/ElasticCableConstantGravity/rho") )
	{
		std::cerr << "Error: Must input density for ElasticCableConstantGravity." << std::endl
				  << "       Please set Physics/ElasticCableConstantGravity/rho." << std::endl;
		libmesh_error();
	}

    int num_gravity =  input.vector_variable_size("Physics/ElasticCableConstantGravity/gravity");
    if (num_gravity != 3)
    {
		std::cerr << "Error: Must input three values for ElasticCableConstantGravity gravity." << std::endl
				  << "       Please set Physics/ElasticCableConstantGravity/gravity." << std::endl;
		libmesh_error();
	}
    for( int i = 0; i < num_gravity; i++ )
	{
    	_gravity(i)=( input("Physics/ElasticCableConstantGravity/gravity", 0.0 , i ) );
	}

    return;
  }

  ElasticCableConstantGravity::~ElasticCableConstantGravity()
  {
    return;
  }

  void ElasticCableConstantGravity::init_variables( libMesh::FEMSystem* system )
  {
    // is_2D = false, is_3D = true
    _disp_vars.init(system,false,true);

    return;
  }


  void ElasticCableConstantGravity::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Tell the system to march temperature forward in time
    system->time_evolving(_disp_vars.u_var());
    system->time_evolving(_disp_vars.v_var());
    system->time_evolving(_disp_vars.w_var());

    return;
  }

  void ElasticCableConstantGravity::init_context( AssemblyContext& context )
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


  void ElasticCableConstantGravity::element_time_derivative( bool compute_jacobian,
                                                                 AssemblyContext& context,
                                                                 CachedValues& /*cache*/ )
  {
    const unsigned int n_u_dofs = context.get_dof_indices(_disp_vars.u_var()).size();

    const std::vector<libMesh::Real> &JxW = context.get_element_fe(_disp_vars.u_var())->get_JxW();

    const std::vector<std::vector<libMesh::Real> >& u_phi = context.get_element_fe(_disp_vars.u_var())->get_phi();

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(_disp_vars.u_var());
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(_disp_vars.v_var());
    libMesh::DenseSubVector<libMesh::Number> &Fw = context.get_elem_residual(_disp_vars.w_var());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    // All shape function gradients are w.r.t. master element coordinates
    const std::vector<std::vector<libMesh::Real> >& dphi_dxi = context.get_element_fe(_disp_vars.u_var())->get_dphidxi();

    const libMesh::DenseSubVector<libMesh::Number>& u_coeffs = context.get_elem_solution( _disp_vars.u_var() );
    const libMesh::DenseSubVector<libMesh::Number>& v_coeffs = context.get_elem_solution( _disp_vars.v_var() );
    const libMesh::DenseSubVector<libMesh::Number>& w_coeffs = context.get_elem_solution( _disp_vars.w_var() );

    //const std::vector<libMesh::RealGradient>& dxdxi  = context.get_element_fe(_disp_vars.u_var())->get_dxyzdxi();

    //Grab the Jacobian matrix
	libMesh::DenseMatrix<libMesh::Number> &K = context.get_elem_jacobian();

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

		libMesh::RealGradient dudxi( grad_u(0), grad_v(0), grad_w(0) );

		libMesh::Real jac = JxW[qp];

		for (unsigned int i=0; i != n_u_dofs; i++)
		{
			Fu(i) += _gravity(0)*_A*_rho*u_phi[i][qp]*jac;

			Fv(i) += _gravity(1)*_A*_rho*u_phi[i][qp]*jac;

			Fw(i) += _gravity(2)*_A*_rho*u_phi[i][qp]*jac;

			/*
			if( compute_jacobian )
			{
				//libmesh_not_implemented();
			}
			*/
		}
	}

    return;
  }
} // end namespace GRINS
