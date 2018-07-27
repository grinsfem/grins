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
#include "grins/axisym_heat_transfer.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/generic_ic_handler.h"
#include "grins/grins_enums.h"
#include "grins/heat_transfer_macros.h"
#include "grins/materials_parsing.h"
#include "grins/variables_parsing.h"
#include "grins/variable_warehouse.h"

// libMesh
#include "libmesh/utility.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace GRINS
{

  template< class Conductivity>
  AxisymmetricHeatTransfer<Conductivity>::AxisymmetricHeatTransfer( const std::string& physics_name,
                                                                    const GetPot& input)
    : Physics(physics_name, input),
      _flow_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<VelocityVariable>(VariablesParsing::velocity_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _press_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PressureFEVariable>(VariablesParsing::press_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _temp_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<PrimitiveTempFEVariables>(VariablesParsing::temp_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _k(input,MaterialsParsing::material_name(input,PhysicsNaming::axisymmetric_heat_transfer()))
  {
    this->read_input_options(input);

    this->_ic_handler = new GenericICHandler( physics_name, input );
  }

  template< class Conductivity>
  void AxisymmetricHeatTransfer<Conductivity>::read_input_options( const GetPot& input )
  {
    MaterialsParsing::read_property( input, "Density", PhysicsNaming::axisymmetric_heat_transfer(),  (*this), this->_rho );
    MaterialsParsing::read_property( input, "SpecificHeat", PhysicsNaming::axisymmetric_heat_transfer(),  (*this), this->_Cp );
  }

  template< class Conductivity>
  void AxisymmetricHeatTransfer<Conductivity>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Tell the system to march temperature forward in time
    system->time_evolving(this->_temp_vars.T(),1);
  }

  template< class Conductivity>
  void AxisymmetricHeatTransfer<Conductivity>::init_context( AssemblyContext& context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    context.get_element_fe(_temp_vars.T())->get_JxW();
    context.get_element_fe(_temp_vars.T())->get_phi();
    context.get_element_fe(_temp_vars.T())->get_dphi();
    context.get_element_fe(_temp_vars.T())->get_xyz();

    context.get_side_fe(_temp_vars.T())->get_JxW();
    context.get_side_fe(_temp_vars.T())->get_phi();
    context.get_side_fe(_temp_vars.T())->get_dphi();
    context.get_side_fe(_temp_vars.T())->get_xyz();
  }

  template< class Conductivity>
  void AxisymmetricHeatTransfer<Conductivity>::element_time_derivative
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.get_dof_indices(_temp_vars.T()).size();
    const unsigned int n_u_dofs = context.get_dof_indices(_flow_vars.u()).size();

    //TODO: check n_T_dofs is same as n_u_dofs, n_v_dofs, n_w_dofs

    // We get some references to cell-specific data that
    // will be used to assemble the linear system.

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(_temp_vars.T())->get_JxW();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(_temp_vars.T())->get_phi();

    // The velocity shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& vel_phi =
      context.get_element_fe(_flow_vars.u())->get_phi();

    // The temperature shape function gradients (in global coords.)
    // at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(_temp_vars.T())->get_dphi();

    // Physical location of the quadrature points
    const std::vector<libMesh::Point>& u_qpoint =
      context.get_element_fe(_flow_vars.u())->get_xyz();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<libMesh::Number> &FT = context.get_elem_residual(_temp_vars.T()); // R_{T}

    libMesh::DenseSubMatrix<libMesh::Number> &KTT = context.get_elem_jacobian(_temp_vars.T(), _temp_vars.T()); // R_{T},{T}

    libMesh::DenseSubMatrix<libMesh::Number> &KTr = context.get_elem_jacobian(_temp_vars.T(), _flow_vars.u()); // R_{T},{r}
    libMesh::DenseSubMatrix<libMesh::Number> &KTz = context.get_elem_jacobian(_temp_vars.T(), _flow_vars.v()); // R_{T},{z}


    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        const libMesh::Number r = u_qpoint[qp](0);

        // Compute the solution & its gradient at the old Newton iterate.
        libMesh::Number u_r, u_z;
        u_r = context.interior_value(_flow_vars.u(), qp);
        u_z = context.interior_value(_flow_vars.v(), qp);

        libMesh::Gradient grad_T;
        grad_T = context.interior_gradient(_temp_vars.T(), qp);

        libMesh::NumberVectorValue U (u_r,u_z);

        libMesh::Number k = this->_k( context, qp );

        // FIXME - once we have T-dependent k, we'll need its
        // derivatives in Jacobians
        // libMesh::Number dk_dT = this->_k.deriv( T );

        // First, an i-loop over the  degrees of freedom.
        for (unsigned int i=0; i != n_T_dofs; i++)
          {
            FT(i) += JxW[qp]*r*
              (-_rho*_Cp*T_phi[i][qp]*(U*grad_T)    // convection term
               -k*(T_gradphi[i][qp]*grad_T) );  // diffusion term

            if (compute_jacobian)
              {
                libmesh_assert (context.get_elem_solution_derivative() == 1.0);

                for (unsigned int j=0; j != n_T_dofs; j++)
                  {
                    // TODO: precompute some terms like:
                    //   _rho*_Cp*T_phi[i][qp]*(vel_phi[j][qp]*T_grad_phi[j][qp])

                    KTT(i,j) += JxW[qp] * context.get_elem_solution_derivative() *r*
                      (-_rho*_Cp*T_phi[i][qp]*(U*T_gradphi[j][qp])  // convection term
                       -k*(T_gradphi[i][qp]*T_gradphi[j][qp])); // diffusion term
                  } // end of the inner dof (j) loop

#if 0
                if( dk_dT != 0.0 )
                  {
                    for (unsigned int j=0; j != n_T_dofs; j++)
                      {
                        // TODO: precompute some terms like:
                        KTT(i,j) -= JxW[qp] * context.get_elem_solution_derivative() *r*( dk_dT*T_phi[j][qp]*T_gradphi[i][qp]*grad_T );
                      }
                  }
#endif

                // Matrix contributions for the Tu, Tv and Tw couplings (n_T_dofs same as n_u_dofs, n_v_dofs and n_w_dofs)
                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    KTr(i,j) += JxW[qp] * context.get_elem_solution_derivative() *r*(-_rho*_Cp*T_phi[i][qp]*(vel_phi[j][qp]*grad_T(0)));
                    KTz(i,j) += JxW[qp] * context.get_elem_solution_derivative() *r*(-_rho*_Cp*T_phi[i][qp]*(vel_phi[j][qp]*grad_T(1)));
                  } // end of the inner dof (j) loop

              } // end - if (compute_jacobian && context.get_elem_solution_derivative())

          } // end of the outer dof (i) loop
      } // end of the quadrature point (qp) loop
  }

  template< class Conductivity>
  void AxisymmetricHeatTransfer<Conductivity>::mass_residual
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // First we get some references to cell-specific data that
    // will be used to assemble the linear system.

    // Element Jacobian * quadrature weights for interior integration
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(_temp_vars.T())->get_JxW();

    // The shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& phi =
      context.get_element_fe(_temp_vars.T())->get_phi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_T_dofs = context.get_dof_indices(_temp_vars.T()).size();

    // Physical location of the quadrature points
    const std::vector<libMesh::Point>& u_qpoint =
      context.get_element_fe(_flow_vars.u())->get_xyz();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<libMesh::Real> &F =
      context.get_elem_residual(_temp_vars.T());

    libMesh::DenseSubMatrix<libMesh::Real> &M =
      context.get_elem_jacobian(_temp_vars.T(), _temp_vars.T());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
        const libMesh::Number r = u_qpoint[qp](0);

        // For the mass residual, we need to be a little careful.
        // The time integrator is handling the time-discretization
        // for us so we need to supply M(u_fixed)*u' for the residual.
        // u_fixed will be given by the fixed_interior_value function
        // while u' will be given by the interior_rate function.
        libMesh::Real T_dot;
        context.interior_rate(_temp_vars.T(), qp, T_dot);

        for (unsigned int i = 0; i != n_T_dofs; ++i)
          {
            F(i) -= JxW[qp]*r*(_rho*_Cp*T_dot*phi[i][qp] );

            if( compute_jacobian )
              {
                for (unsigned int j=0; j != n_T_dofs; j++)
                  {
                    // We're assuming rho, cp are constant w.r.t. T here.
                    M(i,j) -= JxW[qp] * context.get_elem_solution_rate_derivative() *r*_rho*_Cp*phi[j][qp]*phi[i][qp] ;
                  }
              }// End of check on Jacobian

          } // End of element dof loop

      } // End of the quadrature point loop
  }

  template< class Conductivity>
  void AxisymmetricHeatTransfer<Conductivity>::register_parameter
  ( const std::string & param_name,
    libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
    const
  {
    ParameterUser::register_parameter(param_name, param_pointer);
    _k.register_parameter(param_name, param_pointer);
  }


} // namespace GRINS

// Instantiate
INSTANTIATE_HEAT_TRANSFER_SUBCLASS(AxisymmetricHeatTransfer);
