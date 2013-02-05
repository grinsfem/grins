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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "grins/inc_navier_stokes.h"

// libMesh
#include "libmesh/utility.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/fem_context.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  IncompressibleNavierStokesBase::IncompressibleNavierStokesBase(const std::string& physics_name, const GetPot& input )
    : Physics(physics_name, input)
  {
    this->read_input_options(input);

    return;
  }

  IncompressibleNavierStokesBase::~IncompressibleNavierStokesBase()
  {
    return;
  }

  void IncompressibleNavierStokesBase::read_input_options( const GetPot& input )
  {
    // Read FE info
    this->_FE_family =
      libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+incompressible_navier_stokes+"/FE_family", "LAGRANGE") );

    this->_V_order =
      libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+incompressible_navier_stokes+"/V_order", "SECOND") );

    this->_P_order =
      libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+incompressible_navier_stokes+"/P_order", "FIRST") );

    // Read material parameters
    this->_rho = input("Physics/"+incompressible_navier_stokes+"/rho", 1.0);
    this->_mu  = input("Physics/"+incompressible_navier_stokes+"/mu", 1.0);

    // Read variable naming info
    this->_u_var_name = input("Physics/VariableNames/u_velocity", u_var_name_default );
    this->_v_var_name = input("Physics/VariableNames/v_velocity", v_var_name_default );
    this->_w_var_name = input("Physics/VariableNames/w_velocity", w_var_name_default );
    this->_p_var_name = input("Physics/VariableNames/pressure", p_var_name_default );

    return;
  }

  void IncompressibleNavierStokesBase::init_variables( libMesh::FEMSystem* system )
  {
    // Get libMesh to assign an index for each variable
    this->_dim = system->get_mesh().mesh_dimension();

    _u_var = system->add_variable( _u_var_name, this->_V_order, _FE_family);
    _v_var = system->add_variable( _v_var_name, this->_V_order, _FE_family);

    if (_dim == 3)
      _w_var = system->add_variable( _w_var_name, this->_V_order, _FE_family);

    _p_var = system->add_variable( _p_var_name, this->_P_order, _FE_family);

    return;
  }

  void IncompressibleNavierStokesBase::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    const unsigned int dim = system->get_mesh().mesh_dimension();

    // Tell the system to march velocity forward in time, but
    // leave p as a constraint only
    system->time_evolving(_u_var);
    system->time_evolving(_v_var);

    if (dim == 3)
      system->time_evolving(_w_var);

    return;
  }

  void IncompressibleNavierStokesBase::init_context( libMesh::FEMContext &context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    context.element_fe_var[_u_var]->get_JxW();
    context.element_fe_var[_u_var]->get_phi();
    context.element_fe_var[_u_var]->get_dphi();
    context.element_fe_var[_u_var]->get_xyz();

    context.element_fe_var[_p_var]->get_phi();
    context.element_fe_var[_p_var]->get_xyz();

    context.side_fe_var[_u_var]->get_JxW();
    context.side_fe_var[_u_var]->get_phi();
    context.side_fe_var[_u_var]->get_dphi();
    context.side_fe_var[_u_var]->get_xyz();

    return;
  }

} // namespace GRINS
