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
#include "grins/inc_navier_stokes_base.h"

// GRINS
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/utility.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  IncompressibleNavierStokesBase::IncompressibleNavierStokesBase(const std::string& physics_name, const GetPot& input )
    : Physics(physics_name, input),
      _flow_vars(input, incompressible_navier_stokes),
      _rho( input("Physics/"+incompressible_navier_stokes+"/rho", 1.0) ),
      _mu( input("Physics/"+incompressible_navier_stokes+"/mu", 1.0) )
  {
    return;
  }

  IncompressibleNavierStokesBase::~IncompressibleNavierStokesBase()
  {
    return;
  }

  void IncompressibleNavierStokesBase::init_variables( libMesh::FEMSystem* system )
  {
    this->_dim = system->get_mesh().mesh_dimension();

    _flow_vars.init(system);

    return;
  }

  void IncompressibleNavierStokesBase::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    const unsigned int dim = system->get_mesh().mesh_dimension();

    // Tell the system to march velocity forward in time, but
    // leave p as a constraint only
    system->time_evolving(_flow_vars.u_var());
    system->time_evolving(_flow_vars.v_var());

    if (dim == 3)
      system->time_evolving(_flow_vars.w_var());

    return;
  }

  void IncompressibleNavierStokesBase::init_context( AssemblyContext& context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    context.get_element_fe(_flow_vars.u_var())->get_JxW();
    context.get_element_fe(_flow_vars.u_var())->get_phi();
    context.get_element_fe(_flow_vars.u_var())->get_dphi();
    context.get_element_fe(_flow_vars.u_var())->get_xyz();

    context.get_element_fe(_flow_vars.p_var())->get_phi();
    context.get_element_fe(_flow_vars.p_var())->get_xyz();

    context.get_side_fe(_flow_vars.u_var())->get_JxW();
    context.get_side_fe(_flow_vars.u_var())->get_phi();
    context.get_side_fe(_flow_vars.u_var())->get_dphi();
    context.get_side_fe(_flow_vars.u_var())->get_xyz();

    return;
  }

} // namespace GRINS
