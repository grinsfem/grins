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

} // namespace GRINS

// Instantiate
INSTANTIATE_INC_NS_SUBCLASS(IncompressibleNavierStokesBase);
