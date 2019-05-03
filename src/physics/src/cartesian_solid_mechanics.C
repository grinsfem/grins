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
#include "grins/cartesian_solid_mechanics.h"

// GRINS
#include "grins_config.h"
#include "grins/assembly_context.h"
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/quadrature.h"

// C++
#include <sstream>

namespace GRINS
{
  template<unsigned int Dim>
  CartesianSolidMechanics<Dim>::CartesianSolidMechanics( const PhysicsName & physics_name,
                                                         const PhysicsName & core_physics_name,
                                                         const GetPot & input )
    :  SolidMechanicsAbstract<Dim>(physics_name,core_physics_name,input)
  {
    if( this->_disp_vars.dim() != Dim )
      {
        std::stringstream dimss;
        dimss << Dim;

        std::string msg = "ERROR: "+physics_name+" only valid for "+dimss.str()+" dimensions!\n";
        msg += "       Make sure you have three components in your Displacement type variable.\n";
        libmesh_error_msg(msg);
      }
  }

  template<unsigned int Dim>
  void CartesianSolidMechanics<Dim>::init_context( AssemblyContext & context )
  {
    this->get_fe(context)->get_JxW();
    this->get_fe(context)->get_phi();
    this->get_fe(context)->get_dphi();
    this->get_fe(context)->get_xyz();
  }

  template<unsigned int Dim>
  libMesh::Tensor CartesianSolidMechanics<Dim>::form_def_gradient( const libMesh::Gradient & grad_u,
                                                                   const libMesh::Gradient & grad_v ) const
  {
    return libMesh::Tensor( 1.0+grad_u(0), grad_u(1), 0.0,
                            grad_v(0), 1.0+grad_v(1), 0.0,
                            0.0,           0.0,       1.0 );
  }

  template<unsigned int Dim>
  libMesh::Tensor CartesianSolidMechanics<Dim>::form_def_gradient( const libMesh::Gradient & grad_u,
                                                                   const libMesh::Gradient & grad_v,
                                                                   const libMesh::Gradient & grad_w ) const
  {
    return libMesh::Tensor( 1.0+grad_u(0), grad_u(1), grad_u(2),
                            grad_v(0), 1.0+grad_v(1), grad_v(2),
                            grad_w(0), grad_w(1), 1.0+grad_w(2) );
  }

  template<unsigned int Dim>
  void CartesianSolidMechanics<Dim>::mass_residual( bool compute_jacobian, AssemblyContext & context )
  {
    const MultiphysicsSystem & system = context.get_multiphysics_system();

    const int n_u_dofs = context.get_dof_indices(this->_disp_vars.u()).size();

    // We need to extract the corresponding velocity variable.
    // This allows us to use either a FirstOrderUnsteadySolver
    // or a SecondOrderUnsteadySolver. That is, we get back the velocity variable
    // index for FirstOrderUnsteadySolvers or, if it's a SecondOrderUnsteadySolver,
    // this is actually just giving us back the same variable index.

    // If we only wanted to use a SecondOrderUnsteadySolver, then this
    // step would be unnecessary and we would just
    // populate the _u_var, etc. blocks of the residual and Jacobian.
    unsigned int u_dot_var = system.get_second_order_dot_var(this->_disp_vars.u());
    unsigned int v_dot_var = system.get_second_order_dot_var(this->_disp_vars.v());
    unsigned int w_dot_var = libMesh::invalid_uint;
    if(Dim==3)
      w_dot_var = system.get_second_order_dot_var(this->_disp_vars.w());

    // Residuals that we're populating
    libMesh::DenseSubVector<libMesh::Number> & Fu = context.get_elem_residual(u_dot_var);
    libMesh::DenseSubVector<libMesh::Number> & Fv = context.get_elem_residual(v_dot_var);
    libMesh::DenseSubVector<libMesh::Number> * Fw = nullptr;

    if(Dim==3)
      Fw = &context.get_elem_residual(w_dot_var);

    libMesh::DenseSubMatrix<libMesh::Number> & Kuu = context.get_elem_jacobian(u_dot_var, u_dot_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kvv = context.get_elem_jacobian(v_dot_var, v_dot_var);
    libMesh::DenseSubMatrix<libMesh::Number> * Kww = nullptr;

    if(Dim==3)
      Kww = &context.get_elem_jacobian(w_dot_var, w_dot_var);

    const std::vector<libMesh::Real> & JxW = this->get_fe(context)->get_JxW();

    const std::vector<std::vector<libMesh::Real> > & phi = this->get_fe(context)->get_phi();

    int n_qpoints = context.get_element_qrule().n_points();

    for (int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Real u_ddot, v_ddot, w_ddot;
        context.interior_accel(u_dot_var, qp, u_ddot);
        context.interior_accel(v_dot_var, qp, v_ddot);

        if(Dim==3)
          context.interior_accel(w_dot_var, qp, w_ddot);

        for (int i=0; i != n_u_dofs; i++)
          {
            libMesh::Real phi_jac = (this->_rho)*phi[i][qp]*JxW[qp];

            Fu(i) += u_ddot*phi_jac;
            Fv(i) += v_ddot*phi_jac;

            if(Dim==3)
              (*Fw)(i) += w_ddot*phi_jac;

            if (compute_jacobian)
              {
                for (int j=0; j != n_u_dofs; j++)
                  {
                    libMesh::Real jac_term = phi_jac*phi[j][qp];
                    jac_term *= context.get_elem_solution_accel_derivative();

                    Kuu(i,j) += jac_term;
                    Kvv(i,j) += jac_term;

                    if(Dim==3)
                      (*Kww)(i,j) += jac_term;
                  }
              }
          }
      }
  }

  template class CartesianSolidMechanics<2>;
  template class CartesianSolidMechanics<3>;

} // end namespace GRINS
