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
// $Id:$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "boussinesq_buoyancy.h"

void GRINS::BoussinesqBuoyancy::read_input_options( GetPot& input )
{
  // Read variable naming info
  this->_u_var_name = input("Physics/VariableNames/u_velocity", GRINS::u_var_name_default );
  this->_v_var_name = input("Physics/VariableNames/v_velocity", GRINS::v_var_name_default );
  this->_w_var_name = input("Physics/VariableNames/w_velocity", GRINS::w_var_name_default );
  this->_T_var_name = input("Physics/VariableNames/Temperature", GRINS::T_var_name_default );

  _rho_ref = input("Physics/BoussinesqBuoyancy/rho_ref", 1.0);
  _T_ref = input("Physics/BoussinesqBuoyancy/T_ref", 1.0);;
  _beta_T = input("Physics/BoussinesqBuoyancy/beta_T", 1.0);;

  int n_g_comps = input.vector_variable_size("Physics/BoussinesqBuoyancy/g");
  if( n_g_comps != 3 )
    {
      std::cerr << "Error: Must specify 3 gravity components when inputting"
		<< std::endl
		<< "       gravity vector. Found " << n_g_comps
		<< " gravity components."
		<< std::endl;
      libmesh_error();
    }

  _g(1) = input("Physics/BoussinesqBuoyancy/g", 0.0, 0 );
  _g(2) = input("Physics/BoussinesqBuoyancy/g", 0.0, 1 );
  _g(3) = input("Physics/BoussinesqBuoyancy/g", 0.0, 2 );

  return;
}

void GRINS::BoussinesqBuoyancy::init_variables( libMesh::FEMSystem* system )
{
  return;
}

void GRINS::BoussinesqBuoyancy::register_variable_indices(GRINS::VariableMap &global_map)
{
  _u_var = global_map[_u_var_name];
  _v_var = global_map[_v_var_name];

  if (_dim == 3)
    _w_var = global_map[_w_var_name];

  _T_var = global_map[_T_var_name];

  return;
}

bool GRINS::BoussinesqBuoyancy::element_time_derivative( bool request_jacobian,
							 libMesh::DiffContext& context,
							 libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("BoussinesqBuoyancy::element_time_derivative");
#endif
  
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  if (_dim != 3)
    _w_var = _u_var; // for convenience

  // Get residuals
  libMesh::DenseSubVector<Number> &Fu = *c.elem_subresiduals[_u_var]; // R_{u}
  libMesh::DenseSubVector<Number> &Fv = *c.elem_subresiduals[_v_var]; // R_{v}
  libMesh::DenseSubVector<Number> &Fw = *c.elem_subresiduals[_w_var]; // R_{w}

  // Get Jacobians
  libMesh::DenseSubMatrix<Number> &KuT = *c.elem_subjacobians[_u_var][_T_var]; // R_{u},{T}
  libMesh::DenseSubMatrix<Number> &KvT = *c.elem_subjacobians[_v_var][_T_var]; // R_{v},{T}
  libMesh::DenseSubMatrix<Number> &KwT = *c.elem_subjacobians[_w_var][_T_var]; // R_{w},{T}

  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the solution & its gradient at the old Newton iterate.
      libMesh::Number T;
      T = c.interior_value(_T_var, qp);

    }

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("BoussinesqBuoyancy::element_time_derivative");
#endif

  return request_jacobian;
}

void GRINS::BoussinesqBuoyancy::init_context( libMesh::DiffContext &context )
{
  return;
}

bool GRINS::BoussinesqBuoyancy::side_time_derivative( bool request_jacobian,
						      libMesh::DiffContext& context,
						      libMesh::FEMSystem* system )
{
  return request_jacobian;
}

bool GRINS::BoussinesqBuoyancy::element_constraint( bool request_jacobian,
						    libMesh::DiffContext& context,
						    libMesh::FEMSystem* system )
{
  return request_jacobian;
}

bool GRINS::BoussinesqBuoyancy::side_constraint( bool request_jacobian,
						 libMesh::DiffContext& context,
						 libMesh::FEMSystem* system )
{
  return request_jacobian;
}

bool GRINS::BoussinesqBuoyancy::mass_residual( bool request_jacobian,
					       libMesh::DiffContext& context,
					       libMesh::FEMSystem* system )
{
  return request_jacobian;
}

void GRINS::BoussinesqBuoyancy::build_local_variable_map()
{
  return;
}

