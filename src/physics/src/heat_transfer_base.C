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
#include "grins/heat_transfer_base.h"

// GRINS
#include "grins_config.h"

// libMesh
#include "libmesh/utility.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/fem_context.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  HeatTransferBase::HeatTransferBase( const std::string& physics_name, const GetPot& input )
    : Physics(physics_name, input)
  {
    this->read_input_options(input);

    return;
  }

  HeatTransferBase::~HeatTransferBase()
  {
    return;
  }

  void HeatTransferBase::read_input_options( const GetPot& input )
  {
    this->_T_FE_family =
      libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+heat_transfer+"/FE_family", "LAGRANGE") );

    this->_T_order =
      libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+heat_transfer+"/T_order", "SECOND") );

    this->_rho = input("Physics/"+heat_transfer+"/rho", 1.0); //TODO: same as Incompressible NS
    this->_Cp  = input("Physics/"+heat_transfer+"/Cp", 1.0);
    this->_k  = input("Physics/"+heat_transfer+"/k", 1.0);

    this->_T_var_name = input("Physics/VariableNames/Temperature", T_var_name_default );

    // velocity variables. We assume the same element type and order for all velocities.
    this->_V_FE_family = 
      libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+incompressible_navier_stokes+"/FE_family", "LAGRANGE") );

    this->_V_order =
      libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+incompressible_navier_stokes+"/V_order", "SECOND") );

    this->_u_var_name = input("Physics/VariableNames/u_velocity", u_var_name_default );
    this->_v_var_name = input("Physics/VariableNames/v_velocity", v_var_name_default );
    this->_w_var_name = input("Physics/VariableNames/w_velocity", w_var_name_default );

    return;
  }

  void HeatTransferBase::init_variables( libMesh::FEMSystem* system )
  {
    // Get libMesh to assign an index for each variable
    this->_dim = system->get_mesh().mesh_dimension();

    _T_var = system->add_variable( _T_var_name, this->_T_order, _T_FE_family);
 
    // If these are already added, then we just get the index. 
    _u_var = system->add_variable(_u_var_name, _V_order, _V_FE_family );
    _v_var = system->add_variable(_v_var_name, _V_order, _V_FE_family );
    if (_dim == 3)
      _w_var = system->add_variable(_w_var_name, _V_order, _V_FE_family );

    return;
  }

  void HeatTransferBase::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Tell the system to march temperature forward in time
    system->time_evolving(_T_var);

    return;
  }

  void HeatTransferBase::init_context( libMesh::FEMContext& context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    context.element_fe_var[_T_var]->get_JxW();
    context.element_fe_var[_T_var]->get_phi();
    context.element_fe_var[_T_var]->get_dphi();
    context.element_fe_var[_T_var]->get_xyz();

    context.side_fe_var[_T_var]->get_JxW();
    context.side_fe_var[_T_var]->get_phi();
    context.side_fe_var[_T_var]->get_dphi();
    context.side_fe_var[_T_var]->get_xyz();

    return;
  }

} // namespace GRINS
