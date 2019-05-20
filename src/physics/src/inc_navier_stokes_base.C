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
#include "grins/inc_navier_stokes_base.h"

// GRINS
#include "grins/common.h"
#include "grins/assembly_context.h"
#include "grins/physics_naming.h"
#include "grins/inc_nav_stokes_macro.h"
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
  template<class Mu>
  IncompressibleNavierStokesBase<Mu>::IncompressibleNavierStokesBase(const std::string& my_physics_name,
                                                                     const std::string& core_physics_name,
                                                                     const GetPot& input )
    : Physics(my_physics_name, input),
      _flow_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<VelocityVariable>(VariablesParsing::velocity_variable_name(input,core_physics_name,VariablesParsing::PHYSICS))),
      _press_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PressureFEVariable>(VariablesParsing::press_variable_name(input,core_physics_name,VariablesParsing::PHYSICS))),
      _rho(0.0),
      _mu(input,MaterialsParsing::material_name(input,core_physics_name))
  {
    _press_var.set_is_constraint_var(true);

    MaterialsParsing::read_property( input, "Density", core_physics_name,  (*this), this->_rho );

    this->check_var_subdomain_consistency(_flow_vars);
    this->check_var_subdomain_consistency(_press_var);
  }

  template<class Mu>
  libMesh::Real IncompressibleNavierStokesBase<Mu>::get_viscosity_value(AssemblyContext& context, unsigned int qp) const
  {
    return this->_mu(context, qp);
  }

  template<class Mu>
  void IncompressibleNavierStokesBase<Mu>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    const unsigned int dim = system->get_mesh().mesh_dimension();

    // Tell the system to march velocity forward in time, but
    // leave p as a constraint only

    system->time_evolving(_flow_vars.u(),1);

    if (dim > 1)
      system->time_evolving(_flow_vars.v(),1);

    if (dim == 3)
      system->time_evolving(_flow_vars.w(),1);
  }

  template<class Mu>
  void IncompressibleNavierStokesBase<Mu>::init_context( AssemblyContext& context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    context.get_element_fe(_flow_vars.u())->get_JxW();
    context.get_element_fe(_flow_vars.u())->get_phi();
    context.get_element_fe(_flow_vars.u())->get_dphi();
    context.get_element_fe(_flow_vars.u())->get_xyz();

    context.get_element_fe(_press_var.p())->get_phi();
    context.get_element_fe(_press_var.p())->get_xyz();

    context.get_side_fe(_flow_vars.u())->get_JxW();
    context.get_side_fe(_flow_vars.u())->get_phi();
    context.get_side_fe(_flow_vars.u())->get_dphi();
    context.get_side_fe(_flow_vars.u())->get_xyz();
  }

  template<class Mu>
  void IncompressibleNavierStokesBase<Mu>::register_parameter
  ( const std::string & param_name,
    libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
    const
  {
    ParameterUser::register_parameter(param_name, param_pointer);
    _mu.register_parameter(param_name, param_pointer);
  }

  template<class Mu>
  void IncompressibleNavierStokesBase<Mu>::mass_residual_impl
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // Element Jacobian * quadrature weights for interior integration
    // We assume the same for each flow variable
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // The shape functions at interior quadrature points.
    // We assume the same for each flow variable
    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    const std::vector<libMesh::Point>& u_qpoint =
      context.get_element_fe(this->_flow_vars.u())->get_xyz();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<libMesh::Real> & F_u = context.get_elem_residual(this->_flow_vars.u());
    libMesh::DenseSubVector<libMesh::Real> & F_v = context.get_elem_residual(this->_flow_vars.v());
    libMesh::DenseSubVector<libMesh::Real> * F_w = NULL;

    libMesh::DenseSubMatrix<libMesh::Real> & M_uu = context.get_elem_jacobian(this->_flow_vars.u(), this->_flow_vars.u());
    libMesh::DenseSubMatrix<libMesh::Real> & M_vv = context.get_elem_jacobian(this->_flow_vars.v(), this->_flow_vars.v());
    libMesh::DenseSubMatrix<libMesh::Real> * M_ww = NULL;


    if( this->_flow_vars.dim() == 3 )
      {
        F_w  = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}
        M_ww = &context.get_elem_jacobian(this->_flow_vars.w(), this->_flow_vars.w());
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
        // For the mass residual, we need to be a little careful.
        // The time integrator is handling the time-discretization
        // for us so we need to supply M(u_fixed)*u' for the residual.
        // u_fixed will be given by the fixed_interior_value function
        // while u' will be given by the interior_rate function.
        libMesh::Real u_dot, v_dot, w_dot = 0.0;
        context.interior_rate(this->_flow_vars.u(), qp, u_dot);
        context.interior_rate(this->_flow_vars.v(), qp, v_dot);

        if(this->_flow_vars.dim() == 3 )
          context.interior_rate(this->_flow_vars.w(), qp, w_dot);

        const libMesh::Number r = u_qpoint[qp](0);

        libMesh::Real jac = JxW[qp];

        if(Physics::is_axisymmetric())
          jac *= r;

        for (unsigned int i = 0; i != n_u_dofs; ++i)
          {
            F_u(i) -= this->_rho*u_dot*u_phi[i][qp]*jac;
            F_v(i) -= this->_rho*v_dot*u_phi[i][qp]*jac;

            if( this->_flow_vars.dim() == 3 )
              (*F_w)(i) -= this->_rho*w_dot*u_phi[i][qp]*jac;

            if( compute_jacobian )
              {
                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    libMesh::Real value = this->_rho*u_phi[i][qp]*u_phi[j][qp]*jac *
                      context.get_elem_solution_rate_derivative();

                    M_uu(i,j) -= value;
                    M_vv(i,j) -= value;

                    if( this->_flow_vars.dim() == 3)
                      (*M_ww)(i,j) -= value;

                  } // End dof loop
              } // End Jacobian check
          } // End dof loop
      } // End quadrature loop
  }

} // namespace GRINS

// Instantiate
INSTANTIATE_INC_NS_SUBCLASS(IncompressibleNavierStokesBase);
