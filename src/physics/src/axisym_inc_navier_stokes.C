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

#include "axisym_inc_navier_stokes.h"

void GRINS::AxisymmetricIncompNavierStokes::read_input_options( GetPot& input )
{
  // Read FE info
  this->_FE_family =
    libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/AxisymIncompNS/FE_family", "LAGRANGE") );

  this->_V_order =
    libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/AxisymIncompNS/V_order", "SECOND") );

  this->_P_order =
    libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/AxisymIncompNS/P_order", "FIRST") );

  // Read material parameters
  this->_rho = input("Physics/AxisymIncompNS/rho", 1.0);
  this->_mu  = input("Physics/AxisymIncompNS/mu", 1.0);

  // Read variable naming info
  this->_u_r_var_name = input("Physics/VariableNames/r_velocity", GRINS::u_r_var_name_default );
  this->_u_z_var_name = input("Physics/VariableNames/z_velocity", GRINS::u_z_var_name_default );

  this->_p_var_name = input("Physics/VariableNames/pressure", GRINS::p_var_name_default );

  // Read boundary condition info
  /** \todo We ought to be able to put this in the base class somehow so
            that it doesn't have to be rewritten for every physics class.
	    Then, the physics only handles the specifics, e.g. reading
	    in boundary velocities. */
  int num_ids = input.vector_variable_size("Physics/AxisymIncompNS/bc_ids");
  int num_bcs = input.vector_variable_size("Physics/AxisymIncompNS/bc_types");

  if( num_ids != num_bcs )
    {
      std::cerr << "Error: Must specify equal number of boundary ids and boundary conditions"
		<< std::endl;
      libmesh_error();
    }
  
  for( int i = 0; i < num_ids; i++ )
    {
      int bc_id = input("Physics/AxisymIncompNS/bc_ids", -1, i );
      std::string bc_type_in = input("Physics/AxisymIncompNS/bc_types", "NULL", i );

      GRINS::BC_TYPES bc_type = _bound_conds.string_to_enum( bc_type_in );

      _bc_map[bc_id] = bc_type;

      std::stringstream ss;
      ss << bc_id;
      std::string bc_id_string = ss.str();
      
      // Now read in auxillary boundary condition information
      switch(bc_type)
	{
	case(PRESCRIBED_VELOCITY):
	  std::vector<double> vel_in(3,0.0);

	  /* Force the user to specify 3 velocity components regardless of dimension.
	     This should make it easier to keep things correct if we want to have 
	     2D flow not be in the x-y plane. */
	  int n_vel_comps = input.vector_variable_size("Physics/AxisymIncompNS/bound_vel_"+bc_id_string);
	  if( n_vel_comps != 2 )
	    {
	      std::cerr << "Error: Must specify 2 velocity components when inputting"
			<< std::endl
			<< "       prescribed velocities. Found " << n_vel_comps
			<< " velocity components."
			<< std::endl;
	      libmesh_error();
	    }

	  /** \todo Need to unit test this somehow. */
	  vel_in[0] = input("Physics/AxisymIncompNS/bound_vel_"+bc_id_string, 0.0, 0 );
	  vel_in[1] = input("Physics/AxisymIncompNS/bound_vel_"+bc_id_string, 0.0, 1 );

	  _vel_boundary_values[bc_id] = vel_in;

	  break;
	} // End switch(bc_type)
    } // End loop on bc_id

  // Read pressure pinning information
  _pin_pressure = input("Physics/AxisymIncompNS/pin_pressure", true );
  _pin_value = input("Physics/AxisymIncompNS/pin_value", 0.0 );

  unsigned int pin_loc_dim = input.vector_variable_size("Physics/AxisymIncompNS/pin_location");

  // If the user is specifying a pin_location, it had better be 2-dimensional
  if( pin_loc_dim != 2 )
    {
      std::cerr << "Error: pressure pin location must be 2 dimensional"
		<< std::endl;
      libmesh_error();
    }

  _pin_location(0) = input("Physics/AxisymIncompNS/pin_location", 0.0, 0 );
  _pin_location(1) = input("Physics/AxisymIncompNS/pin_location", 0.0, 1 );

  return;
}

void GRINS::AxisymmetricIncompNavierStokes::init_variables( libMesh::FEMSystem* system )
{
  // Get libMesh to assign an index for each variable
  this->_dim = system->get_mesh().mesh_dimension();

  _u_r_var = system->add_variable( _u_r_var_name, this->_V_order, _FE_family);
  _u_z_var = system->add_variable( _u_z_var_name, this->_V_order, _FE_family);

  _p_var = system->add_variable( _p_var_name, this->_P_order, _FE_family);

  // Now build the local map
  this->build_local_variable_map();

  return;
}

void GRINS::AxisymmetricIncompNavierStokes::build_local_variable_map()
{
  _var_map[_u_r_var_name] = _u_r_var;
  _var_map[_u_z_var_name] = _u_z_var;
  
  _var_map[_p_var_name] = _p_var;

  this->_local_variable_map_built = true;

  return;
}

void GRINS::AxisymmetricIncompNavierStokes::set_time_evolving_vars( libMesh::FEMSystem* system )
{
  // Tell the system to march velocity forward in time, but
  // leave p as a constraint only
  system->time_evolving(_u_r_var);
  system->time_evolving(_u_z_var);

  return;
}

void GRINS::AxisymmetricIncompNavierStokes::init_context( libMesh::DiffContext &context )
{
  libMesh::FEMContext &c = libmesh_cast_ref<libMesh::FEMContext&>(context);

  // We should prerequest all the data
  // we will need to build the linear system
  // or evaluate a quantity of interest.
  c.element_fe_var[_u_r_var]->get_JxW();
  c.element_fe_var[_u_r_var]->get_phi();
  c.element_fe_var[_u_r_var]->get_dphi();
  c.element_fe_var[_u_r_var]->get_xyz();

  c.element_fe_var[_p_var]->get_phi();
  c.element_fe_var[_p_var]->get_xyz();

  c.side_fe_var[_u_r_var]->get_JxW();
  c.side_fe_var[_u_r_var]->get_phi();
  c.side_fe_var[_u_r_var]->get_dphi();
  c.side_fe_var[_u_r_var]->get_xyz();

  return;
}

bool GRINS::AxisymmetricIncompNavierStokes::element_time_derivative( bool request_jacobian,
								 libMesh::DiffContext& context,
								 libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("AxisymmetricIncompNavierStokes::element_time_derivative");
#endif

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // The number of local degrees of freedom in each variable.
  const unsigned int n_u_dofs = c.dof_indices_var[_u_r_var].size();
  const unsigned int n_p_dofs = c.dof_indices_var[_p_var].size();

  // Check number of dofs is same for _u_r_var, _u_z_var
  libmesh_assert (n_u_dofs == c.dof_indices_var[_u_z_var].size());

  // We get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[_u_r_var]->get_JxW();

  // The velocity shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& vel_phi =
    c.element_fe_var[_u_r_var]->get_phi();

  // The velocity shape function gradients (in global coords.)
  // at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& vel_gradphi =
    c.element_fe_var[_u_r_var]->get_dphi();

  // The pressure shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& p_phi =
    c.element_fe_var[_p_var]->get_phi();

  // Physical location of the quadrature points
  const std::vector<libMesh::Point>& u_qpoint =
    c.element_fe_var[_u_r_var]->get_xyz();

  // The subvectors and submatrices we need to fill:
  //
  // K_{\alpha \beta} = R_{\alpha},{\beta} = \partial{ R_{\alpha} } / \partial{ {\beta} } (where R denotes residual)
  // e.g., for \alpha = v and \beta = u we get: K{vu} = R_{v},{u}
  // Note that Kp_ur, Kp_uz and Fp comes as constraint.
  libMesh::DenseSubVector<Number> &Fr = *c.elem_subresiduals[_u_r_var]; // R_{r}
  libMesh::DenseSubVector<Number> &Fz = *c.elem_subresiduals[_u_z_var]; // R_{z}

  libMesh::DenseSubMatrix<Number> &Krr = *c.elem_subjacobians[_u_r_var][_u_r_var]; // R_{r},{r}
  libMesh::DenseSubMatrix<Number> &Krz = *c.elem_subjacobians[_u_r_var][_u_z_var]; // R_{r},{z}

  libMesh::DenseSubMatrix<Number> &Kzr = *c.elem_subjacobians[_u_z_var][_u_r_var]; // R_{z},{r}
  libMesh::DenseSubMatrix<Number> &Kzz = *c.elem_subjacobians[_u_z_var][_u_z_var]; // R_{z},{z}

  libMesh::DenseSubMatrix<Number> &Krp = *c.elem_subjacobians[_u_r_var][_p_var]; // R_{r},{p}
  libMesh::DenseSubMatrix<Number> &Kzp = *c.elem_subjacobians[_u_z_var][_p_var]; // R_{z},{p}


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
      libMesh::Number p, u_r, u_z;
      p   = c.interior_value(_p_var, qp);
      u_r = c.interior_value(_u_r_var, qp);
      u_z = c.interior_value(_u_z_var, qp);

      libMesh::Gradient gradu_r, gradu_z;
      gradu_r = c.interior_gradient(_u_r_var, qp);
      gradu_z = c.interior_gradient(_u_z_var, qp);

      libMesh::NumberVectorValue U(u_r, u_z);

      // First, an i-loop over the velocity degrees of freedom.
      // We know that n_u_dofs == n_v_dofs so we can compute contributions
      // for both at the same time.
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fr(i) += JxW[qp]*r*
                   ( - _rho*vel_phi[i][qp]*(U*gradu_r)                  // convection term
                     +   p*( vel_gradphi[i][qp](0) + vel_phi[i][qp]/r ) // pressure term
                     - _mu*( vel_gradphi[i][qp]*gradu_r +               // diffusion term
			     u_r*vel_phi[i][qp]/(r*r) )  );             // diffusion term

          Fz(i) += JxW[qp]*r*
                   ( -_rho*vel_phi[i][qp]*(U*gradu_z)     // convection term
                     + p*vel_gradphi[i][qp](1)            // pressure term
                     - _mu*(vel_gradphi[i][qp]*gradu_z) ); // diffusion term

          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);

              for (unsigned int j=0; j != n_u_dofs; j++)
                {
                  // TODO: precompute some terms like:
                  //   (Uvec*vel_gblgradphivec[j][qp]),
                  //   vel_phi[i][qp]*vel_phi[j][qp],
                  //   (vel_gblgradphivec[i][qp]*vel_gblgradphivec[j][qp])

                  Krr(i,j) += JxW[qp]*r*
                              ( - _rho*vel_phi[i][qp]*(U*vel_gradphi[j][qp])      // convection term
                                - _rho*vel_phi[i][qp]*gradu_r(0)*vel_phi[j][qp]   // convection term
                                - _mu*( vel_gradphi[i][qp]*vel_gradphi[j][qp] +
				        vel_phi[i][qp]*vel_phi[j][qp]/(r*r) )   ); // diffusion term

                  Krz(i,j) += JxW[qp]*r*
		    (-_rho*vel_phi[i][qp]*gradu_r(1)*vel_phi[j][qp]); // convection term

                  Kzz(i,j) += JxW[qp]*r*
                              ( - _rho*vel_phi[i][qp]*(U*vel_gradphi[j][qp])       // convection term
                                - _rho*vel_phi[i][qp]*gradu_z(1)*vel_phi[j][qp]   // convection term 
                                - _mu*(vel_gradphi[i][qp]*vel_gradphi[j][qp])   ); // diffusion term 

                  Kzr(i,j) += JxW[qp]*r*
		    (-_rho*vel_phi[i][qp]*gradu_z(0)*vel_phi[j][qp]); // convection term

                } // end of the inner dof (j) loop

              // Matrix contributions for the up, vp and wp couplings
              for (unsigned int j=0; j != n_p_dofs; j++)
                {
                  Krp(i,j) += JxW[qp]*( r*vel_gradphi[i][qp](0)*p_phi[j][qp] +
				        vel_phi[i][qp]*p_phi[j][qp]  );

                  Kzp(i,j) += JxW[qp]*r*vel_gradphi[i][qp](1)*p_phi[j][qp];

                } // end of the inner dof (j) loop

            } // end - if (request_jacobian && c.elem_solution_derivative)

        } // end of the outer dof (i) loop
    } // end of the quadrature point (qp) loop

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("AxisymmetricIncompNavierStokes::element_time_derivative");
#endif

  return request_jacobian;
}

bool GRINS::AxisymmetricIncompNavierStokes::element_constraint( bool request_jacobian,
								libMesh::DiffContext& context,
								libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("AxisymmetricIncompNavierStokes::element_constraint");
#endif

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // The number of local degrees of freedom in each variable.
  const unsigned int n_u_dofs = c.dof_indices_var[_u_r_var].size();
  const unsigned int n_p_dofs = c.dof_indices_var[_p_var].size();

  // We get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[_u_r_var]->get_JxW();

  const std::vector<std::vector<libMesh::Real> >& vel_phi =
    c.element_fe_var[_u_r_var]->get_phi();

  // The velocity shape function gradients (in global coords.)
  // at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& vel_gradphi =
    c.element_fe_var[_u_r_var]->get_dphi();

  // The pressure shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& p_phi =
    c.element_fe_var[_p_var]->get_phi();

  // Physical location of the quadrature points
  const std::vector<libMesh::Point>& u_qpoint =
    c.element_fe_var[_u_r_var]->get_xyz();

  // The subvectors and submatrices we need to fill:
  libMesh::DenseSubVector<Number> &Fp = *c.elem_subresiduals[_p_var]; // R_{p}

  libMesh::DenseSubMatrix<Number> &Kpr = *c.elem_subjacobians[_p_var][_u_r_var]; // R_{p},{r}
  libMesh::DenseSubMatrix<Number> &Kpz = *c.elem_subjacobians[_p_var][_u_z_var]; // R_{p},{z}

  // Add the constraint given by the continuity equation.
  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      const libMesh::Number r = u_qpoint[qp](0);

      libMesh::Number u_r = c.interior_value(_u_r_var, qp);

      // Compute the velocity gradient at the old Newton iterate.
      libMesh::Gradient gradu_r, gradu_z;
      gradu_r = c.interior_gradient(_u_r_var, qp);
      gradu_z = c.interior_gradient(_u_z_var, qp);
      
      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) += JxW[qp]*(  r*p_phi[i][qp]*(gradu_r(0) + gradu_z(1)) +
			      u_r*p_phi[i][qp] );

          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);

              for (unsigned int j=0; j != n_u_dofs; j++)
                {
                  Kpr(i,j) += JxW[qp]*( r*p_phi[i][qp]*vel_gradphi[j][qp](0) +
					p_phi[i][qp]*vel_phi[j][qp] );

                  Kpz(i,j) += JxW[qp]*r*p_phi[i][qp]*vel_gradphi[j][qp](1);

                } // end of the inner dof (j) loop

            } // end - if (request_jacobian && c.elem_solution_derivative)

        } // end of the outer dof (i) loop
    } // end of the quadrature point (qp) loop


  // Pin p = p_value at p_point
  if( _pin_pressure )
    {
      _bound_conds.pin_value( context, request_jacobian, _p_var,
			      _pin_value, _pin_location );
    }
  

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("AxisymmetricIncompNavierStokes::element_constraint");
#endif

  return request_jacobian;
}


bool GRINS::AxisymmetricIncompNavierStokes::side_time_derivative( bool request_jacobian,
							      libMesh::DiffContext& context,
							      libMesh::FEMSystem* system )
{
  return request_jacobian;
}

bool GRINS::AxisymmetricIncompNavierStokes::side_constraint( bool request_jacobian,
							 libMesh::DiffContext& context,
							 libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("AxisymmetricIncompNavierStokes::side_constraint");
#endif

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  const short int boundary_id =
    system->get_mesh().boundary_info->boundary_id(c.elem, c.side);

  libmesh_assert (boundary_id != libMesh::BoundaryInfo::invalid_id);

  std::map< unsigned int, GRINS::BC_TYPES>::const_iterator 
    bc_map_it = _bc_map.find( boundary_id );

  /* We assume that if you didn't put a boundary id in, then you didn't want to
     set a boundary condition on that boundary. */
  if( bc_map_it != _bc_map.end() )
    {
      switch( bc_map_it->second )
	{
	  // No slip boundary condition
	case GRINS::NO_SLIP:
	  {
	    _bound_conds.apply_dirichlet( context, request_jacobian, _u_r_var, 0.0 );
	    
	    _bound_conds.apply_dirichlet( context, request_jacobian, _u_z_var, 0.0 );
	  }
	  break;

	  // Prescribed constant velocity
	case GRINS::PRESCRIBED_VELOCITY:
	  {
	    _bound_conds.apply_dirichlet( context, request_jacobian, 
					  _u_r_var, _vel_boundary_values[boundary_id][0] );

	    _bound_conds.apply_dirichlet( context, request_jacobian, 
					  _u_z_var, _vel_boundary_values[boundary_id][1] );
	  }
	  break;

	  // Inflow 
	case GRINS::INFLOW:
	  {
	    GRINS::BasePointFuncObj* inflow_func = _bound_funcs[boundary_id];
	    
	    if(!inflow_func)
	      {
		std::cerr << "Error: function not attached for inflow boundary " 
			  << boundary_id << std::endl;
		
		libmesh_error();
	      }

	    std::vector<GRINS::VariableIndex> vars(2);
	    vars[0] = _u_r_var;
	    vars[1] = _u_z_var;
	    
	    std::vector<bool> set_vars(2, false);
	    set_vars[0] = true;
	    set_vars[1] = true;

	    _bound_conds.apply_dirichlet( context, request_jacobian, 
					  vars, set_vars, inflow_func );
	  }
	  break;

	case GRINS::AXISYMMETRIC:
	  {
	    _bound_conds.apply_dirichlet( context, request_jacobian, _u_r_var, 0.0 );
	  }
	  break;
	default:
	  {
	    std::cerr << "Error: Invalid BC type for AxisymmetricIncompNavierStokes."
		      << std::endl;
	    libmesh_error();
	  }

	} // End switch on bc type
    } // End if statement

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("AxisymmetricIncompNavierStokes::side_constraint");
#endif

  return request_jacobian;
}

bool GRINS::AxisymmetricIncompNavierStokes::mass_residual( bool request_jacobian,
						       libMesh::DiffContext& context,
						       libMesh::FEMSystem* system )
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // Element Jacobian * quadrature weights for interior integration
  // We assume the same for each flow variable
  const std::vector<Real> &JxW = 
    c.element_fe_var[_u_r_var]->get_JxW();

  // The shape functions at interior quadrature points.
  // We assume the same for each flow variable
  const std::vector<std::vector<Real> >& u_phi = 
    c.element_fe_var[_u_r_var]->get_phi();

  // Physical location of the quadrature points
  const std::vector<libMesh::Point>& u_qpoint =
    c.element_fe_var[_u_r_var]->get_xyz();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[_u_r_var].size();

  // The subvectors and submatrices we need to fill:
  DenseSubVector<Real> &F_r = *c.elem_subresiduals[_u_r_var];
  DenseSubVector<Real> &F_z = *c.elem_subresiduals[_u_z_var];

  DenseSubMatrix<Real> &M_rr = *c.elem_subjacobians[_u_r_var][_u_r_var];
  DenseSubMatrix<Real> &M_zz = *c.elem_subjacobians[_u_z_var][_u_z_var];

  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp = 0; qp != n_qpoints; ++qp)
    {
      // For the mass residual, we need to be a little careful.
      // The time integrator is handling the time-discretization
      // for us so we need to supply M(u_fixed)*u for the residual.
      // u_fixed will be given by the fixed_interior_* functions
      // while u will be given by the interior_* functions.
      Real u_r_dot = c.interior_value(_u_r_var, qp);
      Real u_z_dot = c.interior_value(_u_z_var, qp);
      
      const libMesh::Number r = u_qpoint[qp](0);

      for (unsigned int i = 0; i != n_u_dofs; ++i)
        {
	  F_r(i) += JxW[qp]*r*_rho*u_r_dot*u_phi[i][qp];
	  F_z(i) += JxW[qp]*r*_rho*u_z_dot*u_phi[i][qp];

	  if( request_jacobian )
              {
		for (unsigned int j=0; j != n_u_dofs; j++)
		  {
		    // Assuming rho is constant w.r.t. u, v, w
		    // and T (if Boussinesq added).
		    Real value = JxW[qp]*r*_rho*u_phi[i][qp]*u_phi[j][qp];

		    M_rr(i,j) += value;
		    M_zz(i,j) += value;

		  } // End dof loop
	      } // End Jacobian check
	} // End dof loop
    } // End quadrature loop

  return request_jacobian;
}
