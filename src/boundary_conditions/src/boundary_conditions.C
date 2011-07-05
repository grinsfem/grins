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
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "boundary_conditions.h"

GRINS::BoundaryConditions::BoundaryConditions( )
{
  return;
}

GRINS::BoundaryConditions::~BoundaryConditions( )
{
  return;
}

void GRINS::BoundaryConditions::apply_dirichlet( libMesh::DiffContext &context,
						 const bool request_jacobian,
						 const VariableIndex var, 
						 const double dbc_value, 
						 const double penalty )
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // The number of local degrees of freedom in u variable.
  const unsigned int n_var_dofs = c.dof_indices_var[var].size();

  // Element Jacobian * quadrature weight for side integration.
  const std::vector<libMesh::Real> &JxW_side = c.side_fe_var[var]->get_JxW();

  // The var shape functions at side quadrature points.
  const std::vector<std::vector<libMesh::Real> >& var_phi_side =
    c.side_fe_var[var]->get_phi();

  // Physical location of the quadrature points on the side.
  const std::vector<libMesh::Point>& var_qpoint = c.side_fe_var[var]->get_xyz();

  libMesh::DenseSubVector<Number> &F_var = *c.elem_subresiduals[var]; // residual
  libMesh::DenseSubMatrix<Number> &K_var = *c.elem_subjacobians[var][var]; // jacobian

  unsigned int n_sidepoints = c.side_qrule->n_points();
  for (unsigned int qp=0; qp != n_sidepoints; qp++)
    {
      // Compute the solution at the old Newton iterate
      libMesh::Number var_value = c.side_value(var, qp);

      for (unsigned int i=0; i != n_var_dofs; i++)
	{
	  F_var(i) += JxW_side[qp] * penalty *
	    ( var_value - dbc_value) * var_phi_side[i][qp];

	  if (request_jacobian)
	    {
	      for (unsigned int j=0; j != n_var_dofs; j++)
		{
		  K_var(i,j) += JxW_side[qp] * penalty *
		    var_phi_side[i][qp] * var_phi_side[j][qp];
		}
	    }
	}
    }

  return;
}

GRINS::BC_TYPES GRINS::BoundaryConditions::string_to_enum( const std::string bc_type )
{
  GRINS::BC_TYPES bc_type_out;

  if( bc_type == "no_slip" )
    bc_type_out = GRINS::NO_SLIP;

  else if( bc_type == "no_flow" )
    bc_type_out = GRINS::NO_FLOW;

  else if( bc_type == "prescribed_vel" )
    bc_type_out = GRINS::PRESCRIBED_VELOCITY;

  else if( bc_type == "inflow" )
    bc_type_out = GRINS::INFLOW;

  else if( bc_type == "outflow" )
    bc_type_out = GRINS::OUTFLOW;

  else if( bc_type == "isothermal_wall" )
    bc_type_out = GRINS::ISOTHERMAL_WALL;

  else if( bc_type == "adiabatic_wall" )
    bc_type_out = GRINS::ADIABATIC_WALL;

  else if( bc_type == "prescribed_heat_flux" )
    bc_type_out = GRINS::PRESCRIBED_HEAT_FLUX;
  
  else
    {
      std::cerr << "Error: Invalid bc_type " << bc_type << std::endl;
      libmesh_error();
    }

  return bc_type_out;
}
