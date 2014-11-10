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
#include "grins/velocity_drag.h"

// GRINS
#include "grins/generic_ic_handler.h"
#include "grins/constant_viscosity.h"
#include "grins/parsed_viscosity.h"
#include "grins/inc_nav_stokes_macro.h"

// libMesh
#include "libmesh/quadrature.h"
#include "libmesh/boundary_info.h"
#include "libmesh/parsed_function.h"
#include "libmesh/zero_function.h"

namespace GRINS
{

  template<class Mu>
  VelocityDrag<Mu>::VelocityDrag( const std::string& physics_name, const GetPot& input )
    : IncompressibleNavierStokesBase<Mu>(physics_name, input)
  {
    this->read_input_options(input);

    return;
  }

  template<class Mu>
  VelocityDrag<Mu>::~VelocityDrag()
  {
    return;
  }

  template<class Mu>
  void VelocityDrag<Mu>::read_input_options( const GetPot& input )
  {
    _exponent = input("Physics/"+velocity_drag+"/exponent", libMesh::Real(2));

    std::string coefficient_function =
      input("Physics/"+velocity_penalty+"/coefficient",
        std::string("0"));

    this->_coefficient.reset
      (new libMesh::ParsedFunction<libMesh::Number>(coefficient_function));

    if (coefficient_function == "0")
      std::cout << "Warning! Zero VelocityDrag specified!" << std::endl;
  }

  template<class Mu>
  void VelocityDrag<Mu>::element_time_derivative( bool compute_jacobian,
					      AssemblyContext& context,
					      CachedValues& /* cache */ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("VelocityDrag::element_time_derivative");
#endif

    // Element Jacobian * quadrature weights for interior integration
    const std::vector<libMesh::Real> &JxW = 
      context.get_element_fe(this->_flow_vars.u_var())->get_JxW();

    // The shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& u_phi = 
      context.get_element_fe(this->_flow_vars.u_var())->get_phi();

    const std::vector<libMesh::Point>& u_qpoint = 
      context.get_element_fe(this->_flow_vars.u_var())->get_xyz();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u_var()).size();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubMatrix<libMesh::Number> &Kuu = context.get_elem_jacobian(this->_flow_vars.u_var(), this->_flow_vars.u_var()); // R_{u},{u}
    libMesh::DenseSubMatrix<libMesh::Number> &Kuv = context.get_elem_jacobian(this->_flow_vars.u_var(), this->_flow_vars.v_var()); // R_{u},{v}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvu = context.get_elem_jacobian(this->_flow_vars.v_var(), this->_flow_vars.u_var()); // R_{v},{u}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvv = context.get_elem_jacobian(this->_flow_vars.v_var(), this->_flow_vars.v_var()); // R_{v},{v}

    libMesh::DenseSubMatrix<libMesh::Number>* Kwu = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kwv = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kww = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kuw = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kvw = NULL;

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u_var()); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v_var()); // R_{v}
    libMesh::DenseSubVector<libMesh::Number>* Fw = NULL;

    if( this->_dim == 3 )
      {
        Kuw = &context.get_elem_jacobian(this->_flow_vars.u_var(), this->_flow_vars.w_var()); // R_{u},{w}
        Kvw = &context.get_elem_jacobian(this->_flow_vars.v_var(), this->_flow_vars.w_var()); // R_{v},{w}

        Kwu = &context.get_elem_jacobian(this->_flow_vars.w_var(), this->_flow_vars.u_var()); // R_{w},{u}
        Kwv = &context.get_elem_jacobian(this->_flow_vars.w_var(), this->_flow_vars.v_var()); // R_{w},{v}
        Kww = &context.get_elem_jacobian(this->_flow_vars.w_var(), this->_flow_vars.w_var()); // R_{w},{w}
        Fw  = &context.get_elem_residual(this->_flow_vars.w_var()); // R_{w}
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // Compute the solution at the old Newton iterate.
        libMesh::Number u, v;
        u = context.interior_value(this->_flow_vars.u_var(), qp);
        v = context.interior_value(this->_flow_vars.v_var(), qp);

        libMesh::NumberVectorValue U(u,v);
        if (this->_dim == 3)
          U(2) = context.interior_value(this->_flow_vars.w_var(), qp); // w

        libMesh::Number Umag = U.size();


        libMesh::Number coeff_val = (*_coefficient)(u_qpoint[qp], context.time);

        libMesh::Number F_coeff = std::pow(Umag, _exponent-1) * -coeff_val;

        libMesh::Number J_coeff = compute_jacobian ? 
                                  std::pow(Umag, _exponent-2) *
                                  -coeff_val * (_exponent-1) :
                                  0;

        libMesh::Real jac = JxW[qp];

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) += jac * F_coeff*U(0)*u_phi[i][qp];
            Fv(i) += jac * F_coeff*U(1)*u_phi[i][qp];

            if( this->_dim == 3 )
              (*Fw)(i) += jac * F_coeff*U(2)*u_phi[i][qp];

	    if( compute_jacobian )
              {
                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    Kuu(i,j) += jac *
                      (F_coeff*u_phi[j][qp] +
                       J_coeff*U(0)*U(0)*u_phi[j][qp]/Umag) *
                      u_phi[i][qp];
                    Kuv(i,j) += jac *
                      (J_coeff*U(0)*U(1)*u_phi[j][qp]/Umag) *
                      u_phi[i][qp];

                    Kvu(i,j) += jac *
                      (J_coeff*U(0)*U(1)*u_phi[j][qp]/Umag) *
                      u_phi[i][qp];
                    Kvv(i,j) += jac *
                      (F_coeff*u_phi[j][qp] +
                       J_coeff*U(1)*U(1)*u_phi[j][qp]/Umag) *
                      u_phi[i][qp];

                    if( this->_dim == 3 )
                      {
                        (*Kuw)(i,j) += jac *
                          (J_coeff*U(0)*U(2)*u_phi[j][qp]/Umag) *
                          u_phi[i][qp];
                        (*Kvw)(i,j) += jac *
                          (J_coeff*U(1)*U(2)*u_phi[j][qp]/Umag) *
                          u_phi[i][qp];

                        (*Kwu)(i,j) += jac *
                          (J_coeff*U(0)*U(2)*u_phi[j][qp]/Umag) *
                          u_phi[i][qp];
                        (*Kwv)(i,j) += jac *
                          (J_coeff*U(1)*U(2)*u_phi[j][qp]/Umag) *
                          u_phi[i][qp];
                        (*Kww)(i,j) += jac *
                          (F_coeff*u_phi[j][qp] +
                           J_coeff*U(2)*U(2)*u_phi[j][qp]/Umag) *
                          u_phi[i][qp];
                      }
                  }
              }
          }
      }


#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("VelocityDrag::element_time_derivative");
#endif

    return;
  }

} // namespace GRINS

// Instantiate
INSTANTIATE_INC_NS_SUBCLASS(VelocityDrag);
