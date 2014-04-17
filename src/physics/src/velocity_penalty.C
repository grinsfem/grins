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


// This class
#include "grins/velocity_penalty.h"

// GRINS
#include "grins/generic_ic_handler.h"

// libMesh
#include "libmesh/quadrature.h"
#include "libmesh/boundary_info.h"
#include "libmesh/parsed_function.h"

namespace GRINS
{

  VelocityPenalty::VelocityPenalty( const std::string& physics_name, const GetPot& input )
    : IncompressibleNavierStokesBase(physics_name, input)
  {
    this->read_input_options(input);

    return;
  }

  VelocityPenalty::~VelocityPenalty()
  {
    return;
  }

  void VelocityPenalty::read_input_options( const GetPot& input )
  {
    std::string penalty_function =
      input("Physics/"+velocity_penalty+"/penalty_function",
        std::string("0"));

    this->normal_vector_function.reset
      (new libMesh::ParsedFunction<Number>(penalty_function));

    std::string base_function =
      input("Physics/"+velocity_penalty+"/base_velocity",
        std::string("0"));

    if (penalty_function == "0" && base_function == "0")
      std::cout << "Warning! Zero VelocityPenalty specified!" << std::endl;

    this->base_velocity_function.reset
      (new libMesh::ParsedFunction<Number>(base_function));

  }

  void VelocityPenalty::element_constraint( bool compute_jacobian,
					    AssemblyContext& context,
					    CachedValues& cache )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("VelocityPenalty::element_constraint");
#endif

    // Element Jacobian * quadrature weights for interior integration
    const std::vector<libMesh::Real> &JxW = 
      context.get_element_fe(this->_flow_vars.u_var())->get_JxW();

    // The shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& u_phi = 
      context.get_element_fe(this->_flow_vars.u_var())->get_phi();

    const std::vector<std::vector<libMesh::Real> >& p_phi = 
      context.get_element_fe(this->_flow_vars.p_var())->get_phi();

    const std::vector<libMesh::Point>& u_qpoint = 
      context.get_element_fe(this->_flow_vars.u_var())->get_xyz();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = context.get_dof_indices(_flow_vars.u_var()).size();
    const unsigned int n_p_dofs = context.get_dof_indices(_flow_vars.p_var()).size();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubMatrix<libMesh::Number> &Kup = context.get_elem_jacobian(_flow_vars.u_var(), _flow_vars.p_var()); // R_{u},{p}
    libMesh::DenseSubMatrix<libMesh::Number> &Kpu = context.get_elem_jacobian(_flow_vars.p_var(), _flow_vars.u_var()); // R_{p},{u}

    libMesh::DenseSubMatrix<libMesh::Number> &Kvp = context.get_elem_jacobian(_flow_vars.v_var(), _flow_vars.p_var()); // R_{v},{p}
    libMesh::DenseSubMatrix<libMesh::Number> &Kpv = context.get_elem_jacobian(_flow_vars.p_var(), _flow_vars.v_var()); // R_{p},{v}

    libMesh::DenseSubMatrix<libMesh::Number>* Kwp = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kpw = NULL;

    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(_flow_vars.p_var()); // R_{p}
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(_flow_vars.u_var()); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(_flow_vars.v_var()); // R_{v}
    libMesh::DenseSubVector<libMesh::Number>* Fw = NULL;

    if( this->_dim == 3 )
      {
        Kpw = &context.get_elem_jacobian(_flow_vars.p_var(), _flow_vars.w_var()); // R_{p},{w}
        Kwp = &context.get_elem_jacobian(_flow_vars.w_var(), _flow_vars.p_var()); // R_{w},{p}
        Fw  = &context.get_elem_residual(_flow_vars.w_var()); // R_{w}
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // Compute the solution & its gradient at the old Newton iterate.
        libMesh::Number p, u, v;
        p = context.interior_value(_flow_vars.p_var(), qp);
        u = context.interior_value(_flow_vars.u_var(), qp);
        v = context.interior_value(_flow_vars.v_var(), qp);

        libMesh::NumberVectorValue U(u,v);
        if (_dim == 3)
          U(2) = context.interior_value(_flow_vars.w_var(), qp); // w

        // Velocity discrepancy (current velocity minus base velocity)
        // normal to constraint plane, scaled by constraint penalty
        // value
        libmesh_assert(normal_vector_function.get());
        libmesh_assert(base_velocity_function.get());

        DenseVector<Number> output_vec(3);

        (*normal_vector_function)(u_qpoint[qp], context.time,
                                  output_vec);

        libMesh::NumberVectorValue U_N(output_vec(0),
                                       output_vec(1),
                                       output_vec(2));

        (*base_velocity_function)(u_qpoint[qp], context.time,
                                  output_vec);

        libMesh::NumberVectorValue U_B(output_vec(0),
                                       output_vec(1),
                                       output_vec(2));

        libMesh::Real jac = JxW[qp];

        for (unsigned int i=0; i != n_p_dofs; i++)
          {
            Fp(i) += jac *
              (((U-U_B)*U_N)*p_phi[i][qp]);

	    if( compute_jacobian )
              {
                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    Kpu(i,j) += jac *
                      ((u_phi[j][qp]*U_N(0))*p_phi[i][qp]);
                    Kpv(i,j) += jac *
                      ((u_phi[j][qp]*U_N(1))*p_phi[i][qp]);

                    if( this->_dim == 3 )
                      {
                        (*Kpw)(i,j) += jac *
                          ((u_phi[j][qp]*U_N(2))*p_phi[i][qp]);
                      }
                  }
              }
          }

        if (true) // make this optional?
          {
            for (unsigned int i=0; i != n_u_dofs; i++)
              {
                Fu(i) += jac *
                  ((u_phi[i][qp]*U_N(0))*p);
                Fv(i) += jac *
                  ((u_phi[i][qp]*U_N(1))*p);
                if( this->_dim == 3 )
                  {
                    (*Fw)(i) += jac *
                      ((u_phi[i][qp]*U_N(2))*p);
                  }

	        if( compute_jacobian )
                  {
                    for (unsigned int j=0; j != n_p_dofs; j++)
                      {
                        Kup(i,j) += jac *
                          ((u_phi[i][qp]*U_N(0))*p_phi[j][qp]);
                        Kvp(i,j) += jac *
                          ((u_phi[i][qp]*U_N(1))*p_phi[j][qp]);

                        if( this->_dim == 3 )
                          {
                            (*Kwp)(i,j) += jac *
                              ((u_phi[i][qp]*U_N(2))*p_phi[j][qp]);
                          }
                      }
                  }
              }
          }
      }


#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("VelocityPenalty::element_constraint");
#endif

    return;
  }

} // namespace GRINS
