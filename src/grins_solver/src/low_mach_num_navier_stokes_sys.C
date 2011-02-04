//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2010 The PECOS Development Team
//
// Please see http://pecos.ices.utexas.edu for more information.
//
// This file is part of GRINS.
//
// GRINS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GRINS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GRINS.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------
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
