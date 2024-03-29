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
#include "grins/velocity_penalty.h"

// GRINS
#include "grins/inc_nav_stokes_macro.h"
#include "grins/spalart_allmaras_viscosity.h"
#include "grins/postprocessed_quantities.h"

// libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{

  template<class Mu>
  VelocityPenalty<Mu>::VelocityPenalty( const std::string& physics_name, const GetPot& input )
    : VelocityPenaltyBase<Mu>(physics_name, input),
    _velocity_penalty_x_index(0),
    _velocity_penalty_y_index(0),
    _velocity_penalty_z_index(0),
    _velocity_penalty_base_x_index(0),
    _velocity_penalty_base_y_index(0),
    _velocity_penalty_base_z_index(0)
  {}

  template<class Mu>
  void VelocityPenalty<Mu>::init_context( AssemblyContext& context )
  {
    context.get_element_fe(this->_flow_vars.u())->get_xyz();
    context.get_element_fe(this->_flow_vars.u())->get_phi();
  }

  template<class Mu>
  void VelocityPenalty<Mu>::register_postprocessing_vars( const GetPot& input,
                                                          PostProcessedQuantities<libMesh::Real>& postprocessing )
  {
    std::string section = "Physics/"+this->_physics_name+"/output_vars";

    std::string vel_penalty = "vel_penalty";
    if (this->_physics_name == "VelocityPenalty2")
      vel_penalty += '2';

    if (this->_physics_name == "VelocityPenalty3")
      vel_penalty += '3';

    if( input.have_variable(section) )
      {
        unsigned int n_vars = input.vector_variable_size(section);

        for( unsigned int v = 0; v < n_vars; v++ )
          {
            std::string name = input(section,"DIE!",v);

            if( name == std::string("velocity_penalty") )
              {
                _velocity_penalty_x_index =
                  postprocessing.register_quantity( vel_penalty+"_x" );

                _velocity_penalty_y_index =
                  postprocessing.register_quantity( vel_penalty+"_y" );

                _velocity_penalty_z_index =
                  postprocessing.register_quantity( vel_penalty+"_z" );
              }
            else if( name == std::string("velocity_penalty_base") )
              {
                _velocity_penalty_base_x_index =
                  postprocessing.register_quantity( vel_penalty+"_base_x" );

                _velocity_penalty_base_y_index =
                  postprocessing.register_quantity( vel_penalty+"_base_y" );

                _velocity_penalty_base_z_index =
                  postprocessing.register_quantity( vel_penalty+"_base_z" );
              }
            else
              {
                std::cerr << "Error: Invalue output_vars value for "+this->_physics_name << std::endl
                          << "       Found " << name << std::endl
                          << "       Acceptable values are: velocity_penalty" << std::endl
                          << "                              velocity_penalty_base" << std::endl;
                libmesh_error();
              }
          }
      }
  }

  template<class Mu>
  void VelocityPenalty<Mu>::element_time_derivative
  ( bool compute_jacobian,
    AssemblyContext & context )
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

    libMesh::DenseSubMatrix<libMesh::Number>* Kwu = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kwv = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kww = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kuw = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kvw = NULL;

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u()); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v()); // R_{v}
    libMesh::DenseSubVector<libMesh::Number>* Fw = NULL;

    if( this->_flow_vars.dim() == 3 )
      {
        Kuw = &context.get_elem_jacobian(this->_flow_vars.u(), this->_flow_vars.w()); // R_{u},{w}
        Kvw = &context.get_elem_jacobian(this->_flow_vars.v(), this->_flow_vars.w()); // R_{v},{w}

        Kwu = &context.get_elem_jacobian(this->_flow_vars.w(), this->_flow_vars.u()); // R_{w},{u}
        Kwv = &context.get_elem_jacobian(this->_flow_vars.w(), this->_flow_vars.v()); // R_{w},{v}
        Kww = &context.get_elem_jacobian(this->_flow_vars.w(), this->_flow_vars.w()); // R_{w},{w}
        Fw  = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // Compute the solution at the old Newton iterate.
        libMesh::Number u, v;
        u = context.interior_value(this->_flow_vars.u(), qp);
        v = context.interior_value(this->_flow_vars.v(), qp);

        libMesh::NumberVectorValue U(u,v);
        if (this->_flow_vars.dim() == 3)
          U(2) = context.interior_value(this->_flow_vars.w(), qp); // w

        libMesh::NumberVectorValue F;
        libMesh::NumberTensorValue dFdU;
        libMesh::NumberTensorValue* dFdU_ptr =
          compute_jacobian ? &dFdU : NULL;
        if (!this->compute_force(u_qpoint[qp], context, U, F, dFdU_ptr))
          continue;

        const libMesh::Real jac = JxW[qp];

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            const libMesh::Number jac_i = jac * u_phi[i][qp];

            Fu(i) += F(0)*jac_i;

            Fv(i) += F(1)*jac_i;
            if( this->_flow_vars.dim() == 3 )
              {
                (*Fw)(i) += F(2)*jac_i;
              }

            if( compute_jacobian )
              {
                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    const libMesh::Number jac_ij = context.get_elem_solution_derivative() * jac_i * u_phi[j][qp];
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
  void VelocityPenalty<Mu>::compute_postprocessed_quantity( unsigned int quantity_index,
                                                            const AssemblyContext& context,
                                                            const libMesh::Point& point,
                                                            libMesh::Real& value )
  {
    libMesh::DenseVector<libMesh::Number> output_vec(3);

    if( quantity_index == this->_velocity_penalty_x_index )
      {
        (*this->normal_vector_function)(context, point, context.time, output_vec);

        value = output_vec(0);
      }
    else if( quantity_index == this->_velocity_penalty_y_index )
      {
        (*this->normal_vector_function)(context, point, context.time, output_vec);

        value = output_vec(1);
      }
    else if( quantity_index == this->_velocity_penalty_z_index )
      {
        (*this->normal_vector_function)(context, point, context.time, output_vec);

        value = output_vec(2);
      }
    else if( quantity_index == this->_velocity_penalty_base_x_index )
      {
        (*this->base_velocity_function)(context, point, context.time, output_vec);

        value = output_vec(0);
      }
    else if( quantity_index == this->_velocity_penalty_base_y_index )
      {
        (*this->base_velocity_function)(context, point, context.time, output_vec);

        value = output_vec(1);
      }
    else if( quantity_index == this->_velocity_penalty_base_z_index )
      {
        (*this->base_velocity_function)(context, point, context.time, output_vec);

        value = output_vec(2);
      }
  }

} // namespace GRINS
// Instantiate
INSTANTIATE_INC_NS_SUBCLASS(VelocityPenalty);
