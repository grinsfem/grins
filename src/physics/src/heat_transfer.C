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
#include "grins/heat_transfer.h"

// GRINS
#include "grins_config.h"
#include "grins/assembly_context.h"
#include "grins/generic_ic_handler.h"
#include "grins/heat_transfer_macros.h"
#include "grins/postprocessed_quantities.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/quadrature.h"
#include "libmesh/boundary_info.h"

namespace GRINS
{

  template<class K>
  HeatTransfer<K>::HeatTransfer( const std::string& physics_name, const GetPot& input )
    : HeatTransferBase<K>(physics_name, PhysicsNaming::heat_transfer(), input),
    _k_index(0)
  {
    this->_ic_handler = new GenericICHandler( physics_name, input );
  }

  template<class K>
  void HeatTransfer<K>::register_postprocessing_vars( const GetPot& input,
                                                      PostProcessedQuantities<libMesh::Real>& postprocessing )
  {
    std::string section = "Physics/"+PhysicsNaming::heat_transfer()+"/output_vars";

    if( input.have_variable(section) )
      {
        unsigned int n_vars = input.vector_variable_size(section);

        for( unsigned int v = 0; v < n_vars; v++ )
          {
            std::string name = input(section,"DIE!",v);

            if( name == std::string("k") )
              {
                this->_k_index = postprocessing.register_quantity( name );
              }
            else
              {
                std::cerr << "Error: Invalid output_vars value for "+PhysicsNaming::heat_transfer() << std::endl
                          << "       Found " << name << std::endl
                          << "       Acceptable values are: k" << std::endl;
                libmesh_error();
              }
          }
      }

    return;
  }

  template<class K>
  void HeatTransfer<K>::element_time_derivative
  ( bool compute_jacobian,
    AssemblyContext & context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();

    //TODO: check n_T_dofs is same as n_u_dofs, n_v_dofs, n_w_dofs

    // We get some references to cell-specific data that
    // will be used to assemble the linear system.

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_temp_vars.T())->get_JxW();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_temp_vars.T())->get_phi();

    // The velocity shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& vel_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    // The temperature shape function gradients (in global coords.)
    // at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(this->_temp_vars.T())->get_dphi();

    const std::vector<libMesh::Point>& u_qpoint =
      context.get_element_fe(this->_flow_vars.u())->get_xyz();

    libMesh::DenseSubMatrix<libMesh::Number> &KTT = context.get_elem_jacobian(this->_temp_vars.T(), this->_temp_vars.T()); // R_{T},{T}

    libMesh::DenseSubMatrix<libMesh::Number> &KTu = context.get_elem_jacobian(this->_temp_vars.T(), this->_flow_vars.u()); // R_{T},{u}
    libMesh::DenseSubMatrix<libMesh::Number> &KTv = context.get_elem_jacobian(this->_temp_vars.T(), this->_flow_vars.v()); // R_{T},{v}
    libMesh::DenseSubMatrix<libMesh::Number>* KTw = NULL;

    libMesh::DenseSubVector<libMesh::Number> &FT = context.get_elem_residual(this->_temp_vars.T()); // R_{T}

    if( this->_flow_vars.dim() == 3 )
      {
        KTw = &context.get_elem_jacobian(this->_temp_vars.T(), this->_flow_vars.w()); // R_{T},{w}
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
        libMesh::Number u, v;
        u = context.interior_value(this->_flow_vars.u(), qp);
        v = context.interior_value(this->_flow_vars.v(), qp);

        libMesh::Gradient grad_T;
        grad_T = context.interior_gradient(this->_temp_vars.T(), qp);

        libMesh::NumberVectorValue U (u,v);
        if (this->_flow_vars.dim() == 3)
          U(2) = context.interior_value(this->_flow_vars.w(), qp);

        const libMesh::Number r = u_qpoint[qp](0);

        libMesh::Real jac = JxW[qp];

        // Compute the conductivity at this qp
        libMesh::Real _k_qp = this->_k(context, qp);

        if(Physics::is_axisymmetric())
          {
            jac *= r;
          }

        // First, an i-loop over the  degrees of freedom.
        for (unsigned int i=0; i != n_T_dofs; i++)
          {
            FT(i) += jac *
              (-this->_rho*this->_Cp*T_phi[i][qp]*(U*grad_T)    // convection term
               -_k_qp*(T_gradphi[i][qp]*grad_T) );  // diffusion term

            if (compute_jacobian)
              {
                for (unsigned int j=0; j != n_T_dofs; j++)
                  {
                    // TODO: precompute some terms like:
                    //   this->_rho*this->_Cp*T_phi[i][qp]*(vel_phi[j][qp]*T_grad_phi[j][qp])

                    KTT(i,j) += jac * context.get_elem_solution_derivative() *
                      (-this->_rho*this->_Cp*T_phi[i][qp]*(U*T_gradphi[j][qp])  // convection term
                       -_k_qp*(T_gradphi[i][qp]*T_gradphi[j][qp])); // diffusion term
                  } // end of the inner dof (j) loop

                // Matrix contributions for the Tu, Tv and Tw couplings (n_T_dofs same as n_u_dofs, n_v_dofs and n_w_dofs)
                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    KTu(i,j) += jac * context.get_elem_solution_derivative()*(-this->_rho*this->_Cp*T_phi[i][qp]*(vel_phi[j][qp]*grad_T(0)));
                    KTv(i,j) += jac * context.get_elem_solution_derivative()*(-this->_rho*this->_Cp*T_phi[i][qp]*(vel_phi[j][qp]*grad_T(1)));
                    if (this->_flow_vars.dim() == 3)
                      (*KTw)(i,j) += jac * context.get_elem_solution_derivative()*(-this->_rho*this->_Cp*T_phi[i][qp]*(vel_phi[j][qp]*grad_T(2)));
                  } // end of the inner dof (j) loop

              } // end - if (compute_jacobian && context.get_elem_solution_derivative())

          } // end of the outer dof (i) loop
      } // end of the quadrature point (qp) loop
  }

  template<class K>
  void HeatTransfer<K>::mass_residual
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // First we get some references to cell-specific data that
    // will be used to assemble the linear system.

    // Element Jacobian * quadrature weights for interior integration
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_temp_vars.T())->get_JxW();

    // The shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& phi =
      context.get_element_fe(this->_temp_vars.T())->get_phi();

    const std::vector<libMesh::Point>& u_qpoint =
      context.get_element_fe(this->_flow_vars.u())->get_xyz();

    // The number of local degrees of freedom in each variable
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<libMesh::Real> &F = context.get_elem_residual(this->_temp_vars.T());

    libMesh::DenseSubMatrix<libMesh::Real> &M = context.get_elem_jacobian(this->_temp_vars.T(), this->_temp_vars.T());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
        // For the mass residual, we need to be a little careful.
        // The time integrator is handling the time-discretization
        // for us so we need to supply M(u_fixed)*u' for the residual.
        // u_fixed will be given by the fixed_interior_value function
        // while u' will be given by the interior_rate function.
        libMesh::Real T_dot;
        context.interior_rate(this->_temp_vars.T(), qp, T_dot);

        const libMesh::Number r = u_qpoint[qp](0);

        libMesh::Real jac = JxW[qp];

        if(Physics::is_axisymmetric())
          {
            jac *= r;
          }

        for (unsigned int i = 0; i != n_T_dofs; ++i)
          {
            F(i) -= this->_rho*this->_Cp*T_dot*phi[i][qp]*jac;

            if( compute_jacobian )
              {
                for (unsigned int j=0; j != n_T_dofs; j++)
                  {
                    // We're assuming rho, cp are constant w.r.t. T here.
                    M(i,j) -= this->_rho*this->_Cp*phi[j][qp]*phi[i][qp]*jac * context.get_elem_solution_rate_derivative();
                  }
              }// End of check on Jacobian

          } // End of element dof loop

      } // End of the quadrature point loop
  }

  template<class K>
  void HeatTransfer<K>::compute_postprocessed_quantity( unsigned int quantity_index,
                                                        const AssemblyContext& context,
                                                        const libMesh::Point& point,
                                                        libMesh::Real& value )
  {
    if( quantity_index == this->_k_index )
      value = this->_k(point, context.get_time());
  }

} // namespace GRINS

// Instantiate
INSTANTIATE_HEAT_TRANSFER_SUBCLASS(HeatTransfer);
