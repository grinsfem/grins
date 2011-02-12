//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - a low Mach number Navier-Stokes Finite-Element Solver
//
// Copyright (C) 2010,2011 The PECOS Development Team
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
// Definitions for the LowMachNumberNavierStokesSystem class.
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "low_mach_num_navier_stokes_sys.h"

void GRINS::LowMachNumberNavierStokesSystem::init_data()
{
  // Setup dummy variable for testing.
  Dummy_var = this->add_variable( "Dummy", libMeshEnums::FIRST );

  // Do the parent's initialization after variables are defined
  FEMSystem::init_data();

  return;
}

void GRINS::LowMachNumberNavierStokesSystem::init_context(
						libMesh::DiffContext &context )
{
  libMesh::FEMContext &c = libmesh_cast_ref<libMesh::FEMContext&>(context);

  // We should prerequest all the data
  // we will need to build the linear system
  // or evaluate a quantity of interest.
  c.element_fe_var[Dummy_var]->get_JxW();
  c.element_fe_var[Dummy_var]->get_phi();
  c.element_fe_var[Dummy_var]->get_dphi();
  c.element_fe_var[Dummy_var]->get_xyz();

  c.side_fe_var[Dummy_var]->get_JxW();
  c.side_fe_var[Dummy_var]->get_phi();
  c.side_fe_var[Dummy_var]->get_dphi();
  c.side_fe_var[Dummy_var]->get_xyz();

  return;
}

bool GRINS::LowMachNumberNavierStokesSystem::element_time_derivative(
						bool request_jacobian,
						libMesh::DiffContext& context )
{
  return request_jacobian;
}

bool GRINS::LowMachNumberNavierStokesSystem::side_time_derivative(
						bool request_jacobian,
						libMesh::DiffContext& context )
{
  return request_jacobian;
}

bool GRINS::LowMachNumberNavierStokesSystem::side_constraint(
						bool request_jacobian,
						libMesh::DiffContext& context )
{
  return request_jacobian;
}

bool GRINS::LowMachNumberNavierStokesSystem::mass_residual(
						bool request_jacobian,
						libMesh::DiffContext& context )
{
  return request_jacobian;
}
