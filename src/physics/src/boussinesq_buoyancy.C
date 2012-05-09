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

#include "boussinesq_buoyancy.h"

GRINS::BoussinesqBuoyancy::BoussinesqBuoyancy( const std::string& physics_name, const GetPot& input )
  : Physics(physics_name)
{
  this->read_input_options(input);
  return;
}

GRINS::BoussinesqBuoyancy::~BoussinesqBuoyancy()
{
  return;
}

void GRINS::BoussinesqBuoyancy::read_input_options( const GetPot& input )
{
  this->_V_FE_family =
    libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+incompressible_navier_stokes+"/FE_family", "LAGRANGE") );

  this->_V_order =
    libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+incompressible_navier_stokes+"/V_order", "SECOND") ); 

  this->_T_FE_family =
    libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+heat_transfer+"/FE_family", "LAGRANGE") );

  this->_T_order =
    libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+heat_transfer+"/T_order", "SECOND") );

  // Read variable naming info
  this->_u_var_name = input("Physics/VariableNames/u_velocity", GRINS::u_var_name_default );
  this->_v_var_name = input("Physics/VariableNames/v_velocity", GRINS::v_var_name_default );
  this->_w_var_name = input("Physics/VariableNames/w_velocity", GRINS::w_var_name_default );
  this->_T_var_name = input("Physics/VariableNames/Temperature", GRINS::T_var_name_default );

  _rho_ref = input("Physics/"+boussinesq_buoyancy+"/rho_ref", 1.0);
  _T_ref = input("Physics/"+boussinesq_buoyancy+"/T_ref", 1.0);;
  _beta_T = input("Physics/"+boussinesq_buoyancy+"/beta_T", 1.0);;

   unsigned int g_dim = input.vector_variable_size("Physics/"+boussinesq_buoyancy+"/g");

  _g(0) = input("Physics/"+boussinesq_buoyancy+"/g", 0.0, 0 );
  _g(1) = input("Physics/"+boussinesq_buoyancy+"/g", 0.0, 1 );
  
  if( g_dim == 3)
    _g(2) = input("Physics/"+boussinesq_buoyancy+"/g", 0.0, 2 );

  return;
}

void GRINS::BoussinesqBuoyancy::init_variables( libMesh::FEMSystem* system )
{
  this->_dim = system->get_mesh().mesh_dimension();

  // If these are already added, then we just get the index. 
  _T_var = system->add_variable( _T_var_name, this->_T_order, _T_FE_family);                                                                                                                                                                 
  _u_var = system->add_variable(_u_var_name, _V_order, _V_FE_family );
  _v_var = system->add_variable(_v_var_name, _V_order, _V_FE_family );
  if (_dim == 3)
    _w_var = system->add_variable(_w_var_name, _V_order, _V_FE_family );

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

  // The number of local degrees of freedom in each variable.
  const unsigned int n_u_dofs = c.dof_indices_var[_u_var].size();
  const unsigned int n_T_dofs = c.dof_indices_var[_T_var].size();

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[_u_var]->get_JxW();

  // The velocity shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& vel_phi =
    c.element_fe_var[_u_var]->get_phi();

  // The temperature shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& T_phi =
    c.element_fe_var[_T_var]->get_phi();

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

      // First, an i-loop over the velocity degrees of freedom.
      // We know that n_u_dofs == n_v_dofs so we can compute contributions
      // for both at the same time.
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
	  Fu(i) += -_rho_ref*_beta_T*(T - _T_ref)*_g(0)*vel_phi[i][qp]*JxW[qp];
	  Fv(i) += -_rho_ref*_beta_T*(T - _T_ref)*_g(1)*vel_phi[i][qp]*JxW[qp];

	  if (_dim == 3)
	    Fw(i) += -_rho_ref*_beta_T*(T - _T_ref)*_g(2)*vel_phi[i][qp]*JxW[qp];

	  if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);
              for (unsigned int j=0; j != n_T_dofs; j++)
		{
		  KuT(i,j) += -_rho_ref*_beta_T*_g(0)*vel_phi[i][qp]*T_phi[j][qp]*JxW[qp];
		  KvT(i,j) += -_rho_ref*_beta_T*_g(1)*vel_phi[i][qp]*T_phi[j][qp]*JxW[qp];

		  if (_dim == 3)
		    KwT(i,j) += -_rho_ref*_beta_T*_g(2)*vel_phi[i][qp]*T_phi[j][qp]*JxW[qp];

		} // End j dof loop
	    } // End request_jacobian check

	} // End i dof loop
    } // End quadrature loop

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("BoussinesqBuoyancy::element_time_derivative");
#endif

  return request_jacobian;
}

void GRINS::BoussinesqBuoyancy::init_context( libMesh::DiffContext& )
{
  return;
}

bool GRINS::BoussinesqBuoyancy::side_time_derivative( bool request_jacobian,
						      libMesh::DiffContext&,
						      libMesh::FEMSystem* )
{
  return request_jacobian;
}

bool GRINS::BoussinesqBuoyancy::element_constraint( bool request_jacobian,
						    libMesh::DiffContext&,
						    libMesh::FEMSystem* )
{
  return request_jacobian;
}

bool GRINS::BoussinesqBuoyancy::side_constraint( bool request_jacobian,
						 libMesh::DiffContext&,
						 libMesh::FEMSystem* )
{
  return request_jacobian;
}

bool GRINS::BoussinesqBuoyancy::mass_residual( bool request_jacobian,
					       libMesh::DiffContext&,
					       libMesh::FEMSystem* )
{
  return request_jacobian;
}


