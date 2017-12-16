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
#include "grins/averaged_turbine.h"

// GRINS
#include "grins/generic_ic_handler.h"
#include "grins/variable_name_defaults.h"
#include "grins/inc_nav_stokes_macro.h"

// libMesh
#include "libmesh/quadrature.h"
#include "libmesh/boundary_info.h"

namespace GRINS
{

  template<class Mu>
  AveragedTurbine<Mu>::AveragedTurbine( const std::string& physics_name, const GetPot& input )
    : AveragedTurbineBase<Mu>(physics_name, input)
  {
    this->_ic_handler = new GenericICHandler( physics_name, input );
  }

  template<class Mu>
  void AveragedTurbine<Mu>::init_context( AssemblyContext& context )
  {
    context.get_element_fe(this->_flow_vars.u())->get_xyz();
    context.get_element_fe(this->_flow_vars.u())->get_phi();
  }


  template<class Mu>
  void AveragedTurbine<Mu>::element_time_derivative
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // The shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    const std::vector<libMesh::Point>& u_qpoint =
      context.get_element_fe(this->_flow_vars.u())->get_xyz();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubMatrix<libMesh::Number> &Kuu = context.get_elem_jacobian(this->_flow_vars.u(), this->_flow_vars.u()); // R_{u},{u}
    libMesh::DenseSubMatrix<libMesh::Number> &Kuv = context.get_elem_jacobian(this->_flow_vars.u(), this->_flow_vars.v()); // R_{u},{v}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvu = context.get_elem_jacobian(this->_flow_vars.v(), this->_flow_vars.u()); // R_{v},{u}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvv = context.get_elem_jacobian(this->_flow_vars.v(), this->_flow_vars.v()); // R_{v},{v}

    libMesh::DenseSubMatrix<libMesh::Number> &Kus =
      context.get_elem_jacobian(this->_flow_vars.u(),
                                this->fan_speed_var()); // R_{u},{s}
    libMesh::DenseSubMatrix<libMesh::Number> &Ksu =
      context.get_elem_jacobian(this->fan_speed_var(),
                                this->_flow_vars.u()); // R_{s},{u}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvs =
      context.get_elem_jacobian(this->_flow_vars.v(),
                                this->fan_speed_var()); // R_{v},{s}
    libMesh::DenseSubMatrix<libMesh::Number> &Ksv =
      context.get_elem_jacobian(this->fan_speed_var(),
                                this->_flow_vars.v()); // R_{s},{v}
    libMesh::DenseSubMatrix<libMesh::Number> &Kss =
      context.get_elem_jacobian(this->fan_speed_var(),
                                this->fan_speed_var()); // R_{s},{s}

    libMesh::DenseSubMatrix<libMesh::Number>* Kwu = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kwv = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kww = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kuw = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kvw = NULL;

    libMesh::DenseSubMatrix<libMesh::Number>* Ksw = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kws = NULL;

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u()); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v()); // R_{v}
    libMesh::DenseSubVector<libMesh::Number>* Fw = NULL;

    libMesh::DenseSubVector<libMesh::Number> &Fs = context.get_elem_residual(this->fan_speed_var()); // R_{s}

    if( this->_flow_vars.dim() == 3 )
      {
        Kuw = &context.get_elem_jacobian(this->_flow_vars.u(), this->_flow_vars.w()); // R_{u},{w}
        Kvw = &context.get_elem_jacobian(this->_flow_vars.v(), this->_flow_vars.w()); // R_{v},{w}

        Kwu = &context.get_elem_jacobian(this->_flow_vars.w(), this->_flow_vars.u()); // R_{w},{u}
        Kwv = &context.get_elem_jacobian(this->_flow_vars.w(), this->_flow_vars.v()); // R_{w},{v}
        Kww = &context.get_elem_jacobian(this->_flow_vars.w(), this->_flow_vars.w()); // R_{w},{w}
        Fw  = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}

        Ksw = &context.get_elem_jacobian(this->fan_speed_var(), this->_flow_vars.w()); // R_{s},{w}
        Kws = &context.get_elem_jacobian(this->_flow_vars.w(), this->fan_speed_var()); // R_{w},{s}

        Fw  = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // Compute the solution at the old Newton iterate.
        libMesh::Number u, v, s;
        u = context.interior_value(this->_flow_vars.u(), qp);
        v = context.interior_value(this->_flow_vars.v(), qp);
        s = context.interior_value(this->fan_speed_var(), qp);

        libMesh::NumberVectorValue U(u,v);
        if (this->_flow_vars.dim() == 3)
          U(2) = context.interior_value(this->_flow_vars.w(), qp); // w

        libMesh::NumberVectorValue U_B_1;
        libMesh::NumberVectorValue F;
        libMesh::NumberTensorValue dFdU;
        libMesh::NumberTensorValue* dFdU_ptr =
          compute_jacobian ? &dFdU : NULL;
        libMesh::NumberVectorValue dFds;
        libMesh::NumberVectorValue* dFds_ptr =
          compute_jacobian ? &dFds : NULL;
        if (!this->compute_force(u_qpoint[qp], context.time, U, s,
                                 U_B_1, F, dFdU_ptr, dFds_ptr))
          continue;

        libMesh::Real jac = JxW[qp];

        // Using this dot product to derive torque *depends* on s=1
        // and U_B_1 corresponding to 1 rad/sec base velocity; this
        // means that the length of U_B_1 is equal to radius.

        // F is the force on the air, so *negative* F is the force on
        // the turbine.
        Fs(0) -= U_B_1 * F * jac;

        if (compute_jacobian)
          {
            Kss(0,0) -= U_B_1 * dFds * jac;

            for (unsigned int j=0; j != n_u_dofs; j++)
              {
                libMesh::Real jac_j = JxW[qp] * u_phi[j][qp];

                for (unsigned int d=0; d != 3; ++d)
                  {
                    Ksu(0,j) -= jac_j * U_B_1(d) * dFdU(d,0);
                    Ksv(0,j) -= jac_j * U_B_1(d) * dFdU(d,1);
                  }

                if (this->_flow_vars.dim() == 3)
                  {
                    for (unsigned int d=0; d != 3; ++d)
                      (*Ksw)(0,j) -= jac_j * U_B_1(d) * dFdU(d,2);
                  }

              } // End j dof loop
          }

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            const libMesh::Number jac_i = jac * u_phi[i][qp];

            Fu(i) += F(0)*jac_i;
            Fv(i) += F(1)*jac_i;

            if( this->_flow_vars.dim() == 3 )
              (*Fw)(i) += F(2)*jac_i;

            if( compute_jacobian )
              {
                Kus(i,0) += dFds(0) * jac_i;
                Kvs(i,0) += dFds(1) * jac_i;
                if( this->_flow_vars.dim() == 3 )
                  (*Kws)(i,0) += dFds(2) * jac_i;

                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    const libMesh::Number jac_ij = jac_i * u_phi[j][qp];
                    Kuu(i,j) += jac_ij * dFdU(0,0);
                    Kuv(i,j) += jac_ij * dFdU(0,1);
                    Kvu(i,j) += jac_ij * dFdU(1,0);
                    Kvv(i,j) += jac_ij * dFdU(1,1);

                    if( this->_flow_vars.dim() == 3 )
                      {
                        (*Kuw)(i,j) += jac_ij * dFdU(0,2);
                        (*Kvw)(i,j) += jac_ij * dFdU(1,2);
                        (*Kwu)(i,j) += jac_ij * dFdU(2,0);
                        (*Kwv)(i,j) += jac_ij * dFdU(2,1);
                        (*Kww)(i,j) += jac_ij * dFdU(2,2);
                      }
                  }
              }
          }
      }
  }



  template<class Mu>
  void AveragedTurbine<Mu>::nonlocal_time_derivative
  (bool compute_jacobian, AssemblyContext & context )
  {
    libMesh::DenseSubMatrix<libMesh::Number> &Kss =
      context.get_elem_jacobian(this->fan_speed_var(), this->fan_speed_var()); // R_{s},{s}

    libMesh::DenseSubVector<libMesh::Number> &Fs =
      context.get_elem_residual(this->fan_speed_var()); // R_{s}

    const std::vector<libMesh::dof_id_type>& dof_indices =
      context.get_dof_indices(this->fan_speed_var());

    const libMesh::Number fan_speed =
      context.get_system().current_solution(dof_indices[0]);

    const libMesh::Number output_torque =
      this->torque_function(libMesh::Point(0), fan_speed);

    Fs(0) += output_torque;

    if (compute_jacobian)
      {
        // FIXME: we should replace this FEM with a hook to the AD fparser stuff
        const libMesh::Number epsilon = 1e-6;
        const libMesh::Number output_torque_deriv =
          (this->torque_function(libMesh::Point(0), fan_speed+epsilon) -
           this->torque_function(libMesh::Point(0), fan_speed-epsilon)) / (2*epsilon);

        Kss(0,0) += output_torque_deriv * context.get_elem_solution_derivative();
      }

    return;
  }



  template<class Mu>
  void AveragedTurbine<Mu>::nonlocal_mass_residual
  ( bool compute_jacobian, AssemblyContext & context )
  {
    libMesh::DenseSubMatrix<libMesh::Number> &Kss =
      context.get_elem_jacobian(this->fan_speed_var(), this->fan_speed_var()); // R_{s},{s}

    libMesh::DenseSubVector<libMesh::Number> &Fs =
      context.get_elem_residual(this->fan_speed_var()); // R_{s}

    const libMesh::DenseSubVector<libMesh::Number> &Us =
      context.get_elem_solution_rate(this->fan_speed_var());

    const libMesh::Number& fan_speed = Us(0);

    Fs(0) -= this->moment_of_inertia * fan_speed;

    if (compute_jacobian)
      {
        Kss(0,0) -= this->moment_of_inertia * context.get_elem_solution_rate_derivative();
      }

    return;
  }


} // namespace GRINS

// Instantiate
INSTANTIATE_INC_NS_SUBCLASS(AveragedTurbine);
