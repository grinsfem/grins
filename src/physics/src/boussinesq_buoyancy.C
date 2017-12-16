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
#include "grins_config.h"
#include "grins/boussinesq_buoyancy.h"

// GRINS
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace GRINS
{

  BoussinesqBuoyancy::BoussinesqBuoyancy( const std::string& physics_name, const GetPot& input )
    : BoussinesqBuoyancyBase(physics_name,input)
  {
    return;
  }

  BoussinesqBuoyancy::~BoussinesqBuoyancy()
  {
    return;
  }

  void BoussinesqBuoyancy::element_time_derivative
  ( bool compute_jacobian,
    AssemblyContext & context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.get_dof_indices(_flow_vars.u()).size();
    const unsigned int n_T_dofs = context.get_dof_indices(_temp_vars.T()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(_flow_vars.u())->get_JxW();

    // The velocity shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& vel_phi =
      context.get_element_fe(_flow_vars.u())->get_phi();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(_temp_vars.T())->get_phi();

    // Get residuals
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(_flow_vars.u()); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(_flow_vars.v()); // R_{v}
    libMesh::DenseSubVector<libMesh::Number>* Fw = NULL;

    // Get Jacobians
    libMesh::DenseSubMatrix<libMesh::Number> &KuT = context.get_elem_jacobian(_flow_vars.u(), _temp_vars.T()); // R_{u},{T}
    libMesh::DenseSubMatrix<libMesh::Number> &KvT = context.get_elem_jacobian(_flow_vars.v(), _temp_vars.T()); // R_{v},{T}
    libMesh::DenseSubMatrix<libMesh::Number>* KwT = NULL;



    if( this->_flow_vars.dim() == 3 )
      {
        Fw  = &context.get_elem_residual(_flow_vars.w()); // R_{w}
        KwT = &context.get_elem_jacobian(_flow_vars.w(), _temp_vars.T()); // R_{w},{T}
      }

    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // Compute the solution & its gradient at the old Newton iterate.
        libMesh::Number T;
        T = context.interior_value(_temp_vars.T(), qp);

        // First, an i-loop over the velocity degrees of freedom.
        // We know that n_u_dofs == n_v_dofs so we can compute contributions
        // for both at the same time.
        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) += -_rho*_beta_T*(T - _T_ref)*_g(0)*vel_phi[i][qp]*JxW[qp];
            Fv(i) += -_rho*_beta_T*(T - _T_ref)*_g(1)*vel_phi[i][qp]*JxW[qp];

            if (this->_flow_vars.dim() == 3)
              (*Fw)(i) += -_rho*_beta_T*(T - _T_ref)*_g(2)*vel_phi[i][qp]*JxW[qp];

            if (compute_jacobian)
              {
                for (unsigned int j=0; j != n_T_dofs; j++)
                  {
                    KuT(i,j) += context.get_elem_solution_derivative() *
                      -_rho*_beta_T*_g(0)*vel_phi[i][qp]*T_phi[j][qp]*JxW[qp];
                    KvT(i,j) += context.get_elem_solution_derivative() *
                      -_rho*_beta_T*_g(1)*vel_phi[i][qp]*T_phi[j][qp]*JxW[qp];

                    if (this->_flow_vars.dim() == 3)
                      (*KwT)(i,j) += context.get_elem_solution_derivative() *
                        -_rho*_beta_T*_g(2)*vel_phi[i][qp]*T_phi[j][qp]*JxW[qp];

                  } // End j dof loop
              } // End compute_jacobian check

          } // End i dof loop
      } // End quadrature loop
  }

} // namespace GRINS
