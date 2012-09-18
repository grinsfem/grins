//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "heat_transfer_base.h"

GRINS::HeatTransferBase::HeatTransferBase( const std::string& physics_name, const GetPot& input )
  : Physics(physics_name, input)
{
  this->read_input_options(input);

  return;
}

GRINS::HeatTransferBase::~HeatTransferBase()
{
  return;
}

void GRINS::HeatTransferBase::read_input_options( const GetPot& input )
{
  this->_T_FE_family =
    libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+heat_transfer+"/FE_family", "LAGRANGE") );

  this->_T_order =
    libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+heat_transfer+"/T_order", "SECOND") );

  this->_rho = input("Physics/"+heat_transfer+"/rho", 1.0); //TODO: same as Incompressible NS
  this->_Cp  = input("Physics/"+heat_transfer+"/Cp", 1.0);
  this->_k  = input("Physics/"+heat_transfer+"/k", 1.0);

  this->_T_var_name = input("Physics/VariableNames/Temperature", GRINS::T_var_name_default );

  // velocity variables. We assume the same element type and order for all velocities.
  this->_V_FE_family = 
    libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+incompressible_navier_stokes+"/FE_family", "LAGRANGE") );

  this->_V_order =
    libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+incompressible_navier_stokes+"/V_order", "SECOND") );

  this->_u_var_name = input("Physics/VariableNames/u_velocity", GRINS::u_var_name_default );
  this->_v_var_name = input("Physics/VariableNames/v_velocity", GRINS::v_var_name_default );
  this->_w_var_name = input("Physics/VariableNames/w_velocity", GRINS::w_var_name_default );

  return;
}

void GRINS::HeatTransferBase::init_variables( libMesh::FEMSystem* system )
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

void GRINS::HeatTransferBase::set_time_evolving_vars( libMesh::FEMSystem* system )
{
  const unsigned int dim = system->get_mesh().mesh_dimension();

  // Tell the system to march temperature forward in time
  system->time_evolving(_T_var);

  return;
}

void GRINS::HeatTransferBase::init_context( libMesh::DiffContext &context )
{
  libMesh::FEMContext &c = libmesh_cast_ref<libMesh::FEMContext&>(context);

  // We should prerequest all the data
  // we will need to build the linear system
  // or evaluate a quantity of interest.
  c.element_fe_var[_T_var]->get_JxW();
  c.element_fe_var[_T_var]->get_phi();
  c.element_fe_var[_T_var]->get_dphi();
  c.element_fe_var[_T_var]->get_xyz();

  c.side_fe_var[_T_var]->get_JxW();
  c.side_fe_var[_T_var]->get_phi();
  c.side_fe_var[_T_var]->get_dphi();
  c.side_fe_var[_T_var]->get_xyz();

  //TODO: _u_var is registered so can we assume things related to _u_var
  //      are available in FEMContext

  return;
}
