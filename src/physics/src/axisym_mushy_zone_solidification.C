//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
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

#include "axisym_mushy_zone_solidification.h"

void GRINS::AxisymmetricMushyZoneSolidification::read_input_options( GetPot& input )
{
  // Read variable naming info
  /* \todo We should move this naming stuff to a singleton class */
  this->_u_r_var_name = input("Physics/VariableNames/r_velocity", GRINS::u_r_var_name_default );
  this->_u_z_var_name = input("Physics/VariableNames/z_velocity", GRINS::u_z_var_name_default );
  this->_T_var_name = input("Physics/VariableNames/Temperature", GRINS::T_var_name_default );

 
  return;
}

void GRINS::AxisymmetricMushyZoneSolidification::init_variables( libMesh::FEMSystem* system )
{
  // No variables to initialize, but we must still "build" the local map.
  // In this case, we have no new variables, so the map will be registered as built.
  this->build_local_variable_map();
  return;
}

void GRINS::AxisymmetricMushyZoneSolidification::register_variable_indices(GRINS::VariableMap &global_map)
{
  _u_r_var = global_map[_u_r_var_name];
  _u_z_var = global_map[_u_z_var_name];

  _T_var = global_map[_T_var_name];

  return;
}

bool GRINS::AxisymmetricMushyZoneSolidification::element_time_derivative( bool request_jacobian,
								     libMesh::DiffContext& context,
								     libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("AxisymmetricMushyZoneSolidification::element_time_derivative");
#endif
  
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // The number of local degrees of freedom in each variable.
  const unsigned int n_u_dofs = c.dof_indices_var[_u_r_var].size();
  const unsigned int n_T_dofs = c.dof_indices_var[_T_var].size();

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[_u_r_var]->get_JxW();

  // The velocity shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& vel_phi =
    c.element_fe_var[_u_r_var]->get_phi();

  // The temperature shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& T_phi =
    c.element_fe_var[_T_var]->get_phi();

  // Physical location of the quadrature points
  const std::vector<libMesh::Point>& u_qpoint =
    c.element_fe_var[_u_r_var]->get_xyz();

  // Get residuals
  libMesh::DenseSubVector<Number> &Fr = *c.elem_subresiduals[_u_r_var]; // R_{r}
  libMesh::DenseSubVector<Number> &Fz = *c.elem_subresiduals[_u_z_var]; // R_{z}

  // Get Jacobians
  libMesh::DenseSubMatrix<Number> &Krr = *c.elem_subjacobians[_u_r_var][_u_r_var]; // R_{r},{r}
  libMesh::DenseSubMatrix<Number> &Kzz = *c.elem_subjacobians[_u_z_var][_u_z_var]; // R_{z},{z}
  libMesh::DenseSubMatrix<Number> &KrT = *c.elem_subjacobians[_u_r_var][_T_var]; // R_{r},{T}
  libMesh::DenseSubMatrix<Number> &KzT = *c.elem_subjacobians[_u_z_var][_T_var]; // R_{z},{T}

  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      const libMesh::Number r = u_qpoint[qp](0);

      // Compute the solution & its gradient at the old Newton iterate.
      const libMesh::Number T = c.interior_value(_T_var, qp);

      // First, an i-loop over the velocity degrees of freedom.
      // We know that n_u_r_dofs == n_u_z_dofs so we can compute contributions
      // for both at the same time.
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
	  const libMesh::Number K_perm = this->compute_K_perm( T );

	  Fr(i) += _mu*( u_r - u_cast_r )/K_perm*vel_phi[i][qp]*r*JxW[qp];
	  Fz(i) += _mu*( u_z - u_cast_z )/K_perm*vel_phi[i][qp]*r*JxW[qp];

	  if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);
              for (unsigned int j=0; j != n_T_dofs; j++)
		{
		  // Convenience
		  const libMesh::Number du = _mu/K_perm*vel_phi[j][qp]*vel_phi[i][qp]*r*JxW[qp];
		  Krr(i,j) += du;
		  Kzz(i,j) += du;

		  const libMesh::Number dK_dT = this->dKperm_dT( T );
		  KrT(i,j) += -_mu*( u_r - u_cast_r )/(K_perm*K_perm)*dK_dT*vel_phi[i][qp]*r*JxW[qp];
		  KzT(i,j) += -_mu*( u_z - u_cast_z )/(K_perm*K_perm)*dK_dT*vel_phi[i][qp]*r*JxW[qp];

		} // End j dof loop
	    } // End request_jacobian check

	} // End i dof loop
    } // End quadrature loop

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("AxisymmetricMushyZoneSolidification::element_time_derivative");
#endif

  return request_jacobian;
}

double GRINS::AxisymmetricMushyZoneSolidification::compute_liquid_phi( const double T )
{
  double phi_l = 0.0;

  if( T >= _T_melt + _delta_T ) return phi_l = 1.0;
  else if( T <= _T_melt - _delta_T ) return phi_l = 0.0;
  else
    {
      const double t = (T - (_T_melt - _delta_T ))/(2.0*_delta_T);
      
      phi_l = t*t*(3.0 - 2.0*t);
    }

  return phi_l;
}

double GRINS::AxisymmetricMushyZoneSolidification:: dphi_dT( const double T )
{
  double dphi_dT = 0;

  if( T >= _T_melt + _delta_T ) return dphi_dT = 0.0;
  else if( T <= _T_melt - _delta_T ) return dphi_dT = 0.0;
  else
    {
      // Minor flop optimization
      const double one_over_two_delta_T = 1.0/(2.0*_delta_T);

      const double t = (T - (_T_melt - _delta_T ))*one_over_two_delta_T;

      const double dt_dT = one_over_two_delta_T;
      const double dphi_dt = (6.0*t - t*t);
      dphi_dT = dphi_dt*dt_dT;
    }

  return dphi_dT;
}

void GRINS::AxisymmetricMushyZoneSolidification::init_context( libMesh::DiffContext &context )
{
  return;
}

bool GRINS::AxisymmetricMushyZoneSolidification::side_time_derivative( bool request_jacobian,
						      libMesh::DiffContext& context,
						      libMesh::FEMSystem* system )
{
  return request_jacobian;
}

bool GRINS::AxisymmetricMushyZoneSolidification::element_constraint( bool request_jacobian,
						    libMesh::DiffContext& context,
						    libMesh::FEMSystem* system )
{
  return request_jacobian;
}

bool GRINS::AxisymmetricMushyZoneSolidification::side_constraint( bool request_jacobian,
						 libMesh::DiffContext& context,
						 libMesh::FEMSystem* system )
{
  return request_jacobian;
}

bool GRINS::AxisymmetricMushyZoneSolidification::mass_residual( bool request_jacobian,
					       libMesh::DiffContext& context,
					       libMesh::FEMSystem* system )
{
  return request_jacobian;
}

void GRINS::AxisymmetricMushyZoneSolidification::build_local_variable_map()
{
  // We only have registered variables in this class so the map is built.
  _local_variable_map_built = true;
  return;
}

