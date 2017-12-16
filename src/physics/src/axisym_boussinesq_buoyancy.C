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
#include "grins/axisym_boussinesq_buoyancy.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/grins_enums.h"
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

  AxisymmetricBoussinesqBuoyancy::AxisymmetricBoussinesqBuoyancy( const std::string& physics_name,
                                                                  const GetPot& input )
    : Physics(physics_name, input),
      _flow_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<VelocityVariable>(VariablesParsing::velocity_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _press_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PressureFEVariable>(VariablesParsing::press_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _temp_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<PrimitiveTempFEVariables>(VariablesParsing::temp_variable_name(input,physics_name,VariablesParsing::PHYSICS)))
  {
    this->read_input_options(input);
  }

  void AxisymmetricBoussinesqBuoyancy::read_input_options( const GetPot& input )
  {
    this->set_parameter
      (_rho, input,
       "Physics/"+PhysicsNaming::axisymmetric_boussinesq_buoyancy()+"/rho_ref", 1.0);

    this->set_parameter
      (_T_ref, input, "Physics/"+PhysicsNaming::axisymmetric_boussinesq_buoyancy()+"/T_ref", 1.0);

    this->set_parameter
      (_beta_T, input, "Physics/"+PhysicsNaming::axisymmetric_boussinesq_buoyancy()+"/beta_T", 1.0);

    _g(0) = input("Physics/"+PhysicsNaming::axisymmetric_boussinesq_buoyancy()+"/g", 0.0, 0 );
    _g(1) = input("Physics/"+PhysicsNaming::axisymmetric_boussinesq_buoyancy()+"/g", 0.0, 1 );
  }

  void AxisymmetricBoussinesqBuoyancy::init_context( AssemblyContext& context )
  {
    context.get_element_fe(_flow_vars.u())->get_JxW();
    context.get_element_fe(_flow_vars.u())->get_phi();
    context.get_element_fe(_flow_vars.u())->get_xyz();
    context.get_element_fe(_temp_vars.T())->get_phi();
  }

  void AxisymmetricBoussinesqBuoyancy::element_time_derivative
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

    // Physical location of the quadrature points
    const std::vector<libMesh::Point>& u_qpoint =
      context.get_element_fe(_flow_vars.u())->get_xyz();

    // Get residuals
    libMesh::DenseSubVector<libMesh::Number> &Fr = context.get_elem_residual(_flow_vars.u()); // R_{r}
    libMesh::DenseSubVector<libMesh::Number> &Fz = context.get_elem_residual(_flow_vars.v()); // R_{z}

    // Get Jacobians
    libMesh::DenseSubMatrix<libMesh::Number> &KrT = context.get_elem_jacobian(_flow_vars.u(), _temp_vars.T()); // R_{r},{T}
    libMesh::DenseSubMatrix<libMesh::Number> &KzT = context.get_elem_jacobian(_flow_vars.v(), _temp_vars.T()); // R_{z},{T}

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
        libMesh::Number T;
        T = context.interior_value(_temp_vars.T(), qp);

        // First, an i-loop over the velocity degrees of freedom.
        // We know that n_u_dofs == n_v_dofs so we can compute contributions
        // for both at the same time.
        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fr(i) += -_rho*_beta_T*(T - _T_ref)*_g(0)*vel_phi[i][qp]*r*JxW[qp];
            Fz(i) += -_rho*_beta_T*(T - _T_ref)*_g(1)*vel_phi[i][qp]*r*JxW[qp];

            if (compute_jacobian && context.get_elem_solution_derivative())
              {
                for (unsigned int j=0; j != n_T_dofs; j++)
                  {
                    const libMesh::Number val =
                      -_rho*_beta_T*vel_phi[i][qp]*T_phi[j][qp]*r*JxW[qp]
                      * context.get_elem_solution_derivative();
                    KrT(i,j) += val*_g(0);
                    KzT(i,j) += val*_g(1);
                  } // End j dof loop
              } // End compute_jacobian check

          } // End i dof loop
      } // End quadrature loop
  }

} // namespace GRINS
