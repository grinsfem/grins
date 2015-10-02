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
    : ElasticCableBase(physics_name,input)
  {
    int num_gravity =  input.vector_variable_size("Physics/"+physics_name+"/gravity");
    if (num_gravity != 3)
      {
        std::cerr << "Error: Must input three values for "+physics_name+" gravity." << std::endl
                  << "       Please set Physics/"+physics_name+"/gravity." << std::endl;
        libmesh_error();
      }
    for( int i = 0; i < num_gravity; i++ )
      {
        _gravity(i)=( input("Physics/"+physics_name+"/gravity", 0.0 , i ) );
      }

    return;
  }

  ElasticCableConstantGravity::~ElasticCableConstantGravity()
  {
    return;
  }


  void ElasticCableConstantGravity::element_time_derivative( bool /*compute_jacobian*/,
                                                             AssemblyContext& context,
                                                             CachedValues& /*cache*/ )
  {
    const unsigned int n_u_dofs = context.get_dof_indices(_disp_vars.u_var()).size();

    const std::vector<libMesh::Real> &JxW = this->get_fe(context)->get_JxW();

    const std::vector<std::vector<libMesh::Real> >& u_phi = this->get_fe(context)->get_phi();

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(_disp_vars.u_var());
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(_disp_vars.v_var());
    libMesh::DenseSubVector<libMesh::Number> &Fw = context.get_elem_residual(_disp_vars.w_var());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Real jac = JxW[qp];

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) += _gravity(0)*_A*_rho*u_phi[i][qp]*jac;

            Fv(i) += _gravity(1)*_A*_rho*u_phi[i][qp]*jac;

            Fw(i) += _gravity(2)*_A*_rho*u_phi[i][qp]*jac;
          }
      }

    return;
  }
} // end namespace GRINS
