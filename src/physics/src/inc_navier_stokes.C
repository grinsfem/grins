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

#include "inc_navier_stokes.h"

void GRINS::IncompressibleNavierStokes::read_input_options( const GetPot& input )
{
  // Read FE info
  this->_FE_family =
    libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/IncompNS/FE_family", "LAGRANGE") );

  this->_V_order =
    libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/IncompNS/V_order", "SECOND") );

  this->_P_order =
    libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/IncompNS/P_order", "FIRST") );

  // Read material parameters
  this->_rho = input("Physics/IncompNS/rho", 1.0);
  this->_mu  = input("Physics/IncompNS/mu", 1.0);

  // Read variable naming info
  this->_u_var_name = input("Physics/VariableNames/u_velocity", GRINS::u_var_name_default );
  this->_v_var_name = input("Physics/VariableNames/v_velocity", GRINS::v_var_name_default );
  this->_w_var_name = input("Physics/VariableNames/w_velocity", GRINS::w_var_name_default );
  this->_p_var_name = input("Physics/VariableNames/pressure", GRINS::p_var_name_default );

  // Read boundary condition info
  /** \todo We ought to be able to put this in the base class somehow so
            that it doesn't have to be rewritten for every physics class.
	    Then, the physics only handles the specifics, e.g. reading
	    in boundary velocities. */
  int num_ids = input.vector_variable_size("Physics/IncompNS/bc_ids");
  int num_bcs = input.vector_variable_size("Physics/IncompNS/bc_types");

  if( num_ids != num_bcs )
    {
      std::cerr << "Error: Must specify equal number of boundary ids and boundary conditions"
		<< std::endl;
      libmesh_error();
    }
  
  for( int i = 0; i < num_ids; i++ )
    {
      int bc_id = input("Physics/IncompNS/bc_ids", -1, i );
      std::string bc_type_in = input("Physics/IncompNS/bc_types", "NULL", i );

      GRINS::BC_TYPES bc_type = _bound_conds.string_to_enum( bc_type_in );

      std::stringstream ss;
      ss << bc_id;
      std::string bc_id_string = ss.str();
      
      // Now read in auxillary boundary condition information
      switch(bc_type)
	{
	case(NO_SLIP):
	  {
	    _dirichlet_bc_map[bc_id] = bc_type;
	  }
	  break;
	case(PRESCRIBED_VELOCITY):
	  {
	    _dirichlet_bc_map[bc_id] = bc_type;
	    std::vector<double> vel_in(3,0.0);
	    
	    /* Force the user to specify 3 velocity components regardless of dimension.
	       This should make it easier to keep things correct if we want to have 
	       2D flow not be in the x-y plane. */
	    int n_vel_comps = input.vector_variable_size("Physics/IncompNS/bound_vel_"+bc_id_string);
	    if( n_vel_comps != 3 )
	      {
		std::cerr << "Error: Must specify 3 velocity components when inputting"
			  << std::endl
			  << "       prescribed velocities. Found " << n_vel_comps
			  << " velocity components."
			  << std::endl;
		libmesh_error();
	      }
	    
	    /** \todo Need to unit test this somehow. */
	    vel_in[0] = input("Physics/IncompNS/bound_vel_"+bc_id_string, 0.0, 0 );
	    vel_in[1] = input("Physics/IncompNS/bound_vel_"+bc_id_string, 0.0, 1 );
	    vel_in[2] = input("Physics/IncompNS/bound_vel_"+bc_id_string, 0.0, 2 );
	    
	    _vel_boundary_values[bc_id] = vel_in;
	  }
	  break;
	case(INFLOW):
	  {
	    _dirichlet_bc_map[bc_id] = bc_type;
	  }
	  break;
	default:
	  {
	    std::cerr << "Error: Invalid boundary condition type for IncompressibleNavierStokes."
		      << std::endl;
	    libmesh_error();
	  }
	} // End switch(bc_type)
    } // End loop on bc_id

  // Read pressure pinning information
  _pin_pressure = input("Physics/IncompNS/pin_pressure", true );
  _pin_value = input("Physics/IncompNS/pin_value", 0.0 );

  unsigned int pin_loc_dim = input.vector_variable_size("Physics/IncompNS/pin_location");

  // If the user is specifying a pin_location, it had better be at least 2-dimensional
  if( pin_loc_dim > 0 && pin_loc_dim < 2 )
    {
      std::cerr << "Error: pressure pin location must be at least 2 dimensional"
		<< std::endl;
      libmesh_error();
    }

  _pin_location(0) = input("Physics/IncompNS/pin_location", 0.0, 0 );
  _pin_location(1) = input("Physics/IncompNS/pin_location", 0.0, 1 );

  if( pin_loc_dim == 3 ) 
    _pin_location(2) = input("Physics/IncompNS/pin_location", 0.0, 2 );

  return;
}

void GRINS::IncompressibleNavierStokes::init_variables( libMesh::FEMSystem* system )
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

void GRINS::IncompressibleNavierStokes::set_time_evolving_vars( libMesh::FEMSystem* system )
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

void GRINS::IncompressibleNavierStokes::init_context( libMesh::DiffContext &context )
{
  libMesh::FEMContext &c = libmesh_cast_ref<libMesh::FEMContext&>(context);

  // We should prerequest all the data
  // we will need to build the linear system
  // or evaluate a quantity of interest.
  c.element_fe_var[_u_var]->get_JxW();
  c.element_fe_var[_u_var]->get_phi();
  c.element_fe_var[_u_var]->get_dphi();
  c.element_fe_var[_u_var]->get_xyz();

  c.element_fe_var[_p_var]->get_phi();
  c.element_fe_var[_p_var]->get_xyz();

  c.side_fe_var[_u_var]->get_JxW();
  c.side_fe_var[_u_var]->get_phi();
  c.side_fe_var[_u_var]->get_dphi();
  c.side_fe_var[_u_var]->get_xyz();

  return;
}

bool GRINS::IncompressibleNavierStokes::element_time_derivative( bool request_jacobian,
								 libMesh::DiffContext& context,
								 libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("IncompressibleNavierStokes::element_time_derivative");
#endif

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // The number of local degrees of freedom in each variable.
  const unsigned int n_u_dofs = c.dof_indices_var[_u_var].size();
  const unsigned int n_p_dofs = c.dof_indices_var[_p_var].size();

  // Check number of dofs is same for _u_var, v_var and w_var.
  libmesh_assert (n_u_dofs == c.dof_indices_var[_v_var].size());
  if (_dim == 3)
    libmesh_assert (n_u_dofs == c.dof_indices_var[_w_var].size());

  // We get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[_u_var]->get_JxW();

  // The velocity shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& vel_phi =
    c.element_fe_var[_u_var]->get_phi();

  // The velocity shape function gradients (in global coords.)
  // at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& vel_gblgradphivec =
    c.element_fe_var[_u_var]->get_dphi();

  // The pressure shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& p_phi =
    c.element_fe_var[_p_var]->get_phi();

  // Physical location of the quadrature points
  const std::vector<libMesh::Point>& u_qpoint =
    c.element_fe_var[_u_var]->get_xyz();

  // The subvectors and submatrices we need to fill:
  //
  // K_{\alpha \beta} = R_{\alpha},{\beta} = \partial{ R_{\alpha} } / \partial{ {\beta} } (where R denotes residual)
  // e.g., for \alpha = v and \beta = u we get: K{vu} = R_{v},{u}
  // Note that Kpu, Kpv, Kpw and Fp comes as constraint.
  //
  if (_dim != 3)
    _w_var = _u_var; // for convenience

  libMesh::DenseSubMatrix<Number> &Kuu = *c.elem_subjacobians[_u_var][_u_var]; // R_{u},{u}
  libMesh::DenseSubMatrix<Number> &Kuv = *c.elem_subjacobians[_u_var][_v_var]; // R_{u},{v}
  libMesh::DenseSubMatrix<Number> &Kuw = *c.elem_subjacobians[_u_var][_w_var]; // R_{u},{w}

  libMesh::DenseSubMatrix<Number> &Kvu = *c.elem_subjacobians[_v_var][_u_var]; // R_{v},{u}
  libMesh::DenseSubMatrix<Number> &Kvv = *c.elem_subjacobians[_v_var][_v_var]; // R_{v},{v}
  libMesh::DenseSubMatrix<Number> &Kvw = *c.elem_subjacobians[_v_var][_w_var]; // R_{v},{w}

  libMesh::DenseSubMatrix<Number> &Kwu = *c.elem_subjacobians[_w_var][_u_var]; // R_{w},{u}
  libMesh::DenseSubMatrix<Number> &Kwv = *c.elem_subjacobians[_w_var][_v_var]; // R_{w},{v}
  libMesh::DenseSubMatrix<Number> &Kww = *c.elem_subjacobians[_w_var][_w_var]; // R_{w},{w}

  libMesh::DenseSubMatrix<Number> &Kup = *c.elem_subjacobians[_u_var][_p_var]; // R_{u},{p}
  libMesh::DenseSubMatrix<Number> &Kvp = *c.elem_subjacobians[_v_var][_p_var]; // R_{v},{p}
  libMesh::DenseSubMatrix<Number> &Kwp = *c.elem_subjacobians[_w_var][_p_var]; // R_{w},{p}

  libMesh::DenseSubVector<Number> &Fu = *c.elem_subresiduals[_u_var]; // R_{u}
  libMesh::DenseSubVector<Number> &Fv = *c.elem_subresiduals[_v_var]; // R_{v}
  libMesh::DenseSubVector<Number> &Fw = *c.elem_subresiduals[_w_var]; // R_{w}

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
      libMesh::Number p, u, v, w;
      p = c.interior_value(_p_var, qp);
      u = c.interior_value(_u_var, qp);
      v = c.interior_value(_v_var, qp);
      if (_dim == 3)
        w = c.interior_value(_w_var, qp);

      libMesh::Gradient graduvec, gradvvec, gradwvec;
      graduvec = c.interior_gradient(_u_var, qp);
      gradvvec = c.interior_gradient(_v_var, qp);
      if (_dim == 3)
        gradwvec = c.interior_gradient(_w_var, qp);

      libMesh::NumberVectorValue Uvec (u,v);
      if (_dim == 3)
        Uvec(2) = w;

      const libMesh::Number  graduvec_x = graduvec(0);
      const libMesh::Number  graduvec_y = graduvec(1);
      const libMesh::Number  graduvec_z = (_dim == 3)?graduvec(2):0;
      const libMesh::Number  gradvvec_x = gradvvec(0);
      const libMesh::Number  gradvvec_y = gradvvec(1);
      const libMesh::Number  gradvvec_z = (_dim == 3)?gradvvec(2):0;
      const libMesh::Number  gradwvec_x = (_dim == 3)?gradwvec(0):0;
      const libMesh::Number  gradwvec_y = (_dim == 3)?gradwvec(1):0;
      const libMesh::Number  gradwvec_z = (_dim == 3)?gradwvec(2):0;

      // First, an i-loop over the velocity degrees of freedom.
      // We know that n_u_dofs == n_v_dofs so we can compute contributions
      // for both at the same time.
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fu(i) += JxW[qp] *
                   (-_rho*vel_phi[i][qp]*(Uvec*graduvec)        // convection term
                    +p*vel_gblgradphivec[i][qp](0)              // pressure term
                    -_mu*(vel_gblgradphivec[i][qp]*graduvec) ); // diffusion term

          Fv(i) += JxW[qp] *
                   (-_rho*vel_phi[i][qp]*(Uvec*gradvvec)        // convection term
                    +p*vel_gblgradphivec[i][qp](1)              // pressure term
                    -_mu*(vel_gblgradphivec[i][qp]*gradvvec) ); // diffusion term
          if (_dim == 3)
            {
              Fw(i) += JxW[qp] *
                       (-_rho*vel_phi[i][qp]*(Uvec*gradwvec)        // convection term
                        +p*vel_gblgradphivec[i][qp](2)              // pressure term
                        -_mu*(vel_gblgradphivec[i][qp]*gradwvec) ); // diffusion term
            }

          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);

              for (unsigned int j=0; j != n_u_dofs; j++)
                {
                  // TODO: precompute some terms like:
                  //   (Uvec*vel_gblgradphivec[j][qp]),
                  //   vel_phi[i][qp]*vel_phi[j][qp],
                  //   (vel_gblgradphivec[i][qp]*vel_gblgradphivec[j][qp])

                  Kuu(i,j) += JxW[qp] *
                              (-_rho*vel_phi[i][qp]*(Uvec*vel_gblgradphivec[j][qp])       // convection term
                               -_rho*vel_phi[i][qp]*graduvec_x*vel_phi[j][qp]             // convection term
                               -_mu*(vel_gblgradphivec[i][qp]*vel_gblgradphivec[j][qp])); // diffusion term
                  Kuv(i,j) += JxW[qp] *
                              (-_rho*vel_phi[i][qp]*graduvec_y*vel_phi[j][qp]);           // convection term

                  Kvv(i,j) += JxW[qp] *
                              (-_rho*vel_phi[i][qp]*(Uvec*vel_gblgradphivec[j][qp])       // convection term
                               -_rho*vel_phi[i][qp]*gradvvec_y*vel_phi[j][qp]             // convection term
                               -_mu*(vel_gblgradphivec[i][qp]*vel_gblgradphivec[j][qp])); // diffusion term
                  Kvu(i,j) += JxW[qp] *
                              (-_rho*vel_phi[i][qp]*gradvvec_x*vel_phi[j][qp]);           // convection term

                  if (_dim == 3)
                    {
                      Kuw(i,j) += JxW[qp] *
                                  (-_rho*vel_phi[i][qp]*graduvec_z*vel_phi[j][qp]);           // convection term

                      Kvw(i,j) += JxW[qp] *
                                  (-_rho*vel_phi[i][qp]*gradvvec_z*vel_phi[j][qp]);           // convection term

                      Kww(i,j) += JxW[qp] *
                                  (-_rho*vel_phi[i][qp]*(Uvec*vel_gblgradphivec[j][qp])       // convection term
                                   -_rho*vel_phi[i][qp]*gradwvec_z*vel_phi[j][qp]             // convection term
                                   -_mu*(vel_gblgradphivec[i][qp]*vel_gblgradphivec[j][qp])); // diffusion term
                      Kwu(i,j) += JxW[qp] *
                                  (-_rho*vel_phi[i][qp]*gradwvec_x*vel_phi[j][qp]);           // convection term
                      Kwv(i,j) += JxW[qp] *
                                  (-_rho*vel_phi[i][qp]*gradwvec_y*vel_phi[j][qp]);           // convection term
                    }
                } // end of the inner dof (j) loop

              // Matrix contributions for the up, vp and wp couplings
              for (unsigned int j=0; j != n_p_dofs; j++)
                {
                  Kup(i,j) += JxW[qp]*vel_gblgradphivec[i][qp](0)*p_phi[j][qp];
                  Kvp(i,j) += JxW[qp]*vel_gblgradphivec[i][qp](1)*p_phi[j][qp];
                  if (_dim == 3)
                    Kwp(i,j) += JxW[qp]*vel_gblgradphivec[i][qp](2)*p_phi[j][qp];
                } // end of the inner dof (j) loop

            } // end - if (request_jacobian && c.elem_solution_derivative)

        } // end of the outer dof (i) loop
    } // end of the quadrature point (qp) loop

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("IncompressibleNavierStokes::element_time_derivative");
#endif

  return request_jacobian;
}

bool GRINS::IncompressibleNavierStokes::element_constraint( bool request_jacobian,
							    libMesh::DiffContext& context,
							    libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("IncompressibleNavierStokes::element_constraint");
#endif

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // The number of local degrees of freedom in each variable.
  const unsigned int n_u_dofs = c.dof_indices_var[_u_var].size();
  const unsigned int n_p_dofs = c.dof_indices_var[_p_var].size();

  // We get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[_u_var]->get_JxW();

  // The velocity shape function gradients (in global coords.)
  // at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& vel_gblgradphivec =
    c.element_fe_var[_u_var]->get_dphi();

  // The pressure shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& p_phi =
    c.element_fe_var[_p_var]->get_phi();

  // The subvectors and submatrices we need to fill:
  //
  // Kpu, Kpv, Kpw, Fp
  //
  if (_dim != 3)
    _w_var = _u_var; // for convenience

  libMesh::DenseSubMatrix<Number> &Kpu = *c.elem_subjacobians[_p_var][_u_var]; // R_{p},{u}
  libMesh::DenseSubMatrix<Number> &Kpv = *c.elem_subjacobians[_p_var][_v_var]; // R_{p},{v}
  libMesh::DenseSubMatrix<Number> &Kpw = *c.elem_subjacobians[_p_var][_w_var]; // R_{p},{w}

  libMesh::DenseSubVector<Number> &Fp = *c.elem_subresiduals[_p_var]; // R_{p}

  // Add the constraint given by the continuity equation.
  unsigned int n_qpoints = c.element_qrule->n_points();
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the velocity gradient at the old Newton iterate.
      libMesh::Gradient graduvec, gradvvec, gradwvec;
      graduvec = c.interior_gradient(_u_var, qp);
      gradvvec = c.interior_gradient(_v_var, qp);
      if (_dim == 3)
        gradwvec = c.interior_gradient(_w_var, qp);

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) += JxW[qp] * p_phi[i][qp] *
                   (graduvec(0) + gradvvec(1));
          if (_dim == 3)
            Fp(i) += JxW[qp] * p_phi[i][qp] *
                     (gradwvec(2));

          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);

              for (unsigned int j=0; j != n_u_dofs; j++)
                {
                  Kpu(i,j) += JxW[qp]*p_phi[i][qp]*vel_gblgradphivec[j][qp](0);
                  Kpv(i,j) += JxW[qp]*p_phi[i][qp]*vel_gblgradphivec[j][qp](1);
                  if (_dim == 3)
                    Kpw(i,j) += JxW[qp]*p_phi[i][qp]*vel_gblgradphivec[j][qp](2);
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
  this->_timer->EndTimer("IncompressibleNavierStokes::element_constraint");
#endif

  return request_jacobian;
}


bool GRINS::IncompressibleNavierStokes::side_time_derivative( bool request_jacobian,
							      libMesh::DiffContext& context,
							      libMesh::FEMSystem* system )
{
  return request_jacobian;
}

bool GRINS::IncompressibleNavierStokes::side_constraint( bool request_jacobian,
							 libMesh::DiffContext& context,
							 libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("IncompressibleNavierStokes::side_constraint");
#endif

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // for convenience
  if (_dim != 3)
    _w_var = _u_var;

  const GRINS::BoundaryID boundary_id =
    system->get_mesh().boundary_info->boundary_id(c.elem, c.side);
  libmesh_assert (boundary_id != libMesh::BoundaryInfo::invalid_id);

  std::map< GRINS::BoundaryID, GRINS::BC_TYPES>::const_iterator 
    bc_map_it = _dirichlet_bc_map.find( boundary_id );

  /* We assume that if you didn't put a boundary id in, then you didn't want to
     set a boundary condition on that boundary. */
  if( bc_map_it != _dirichlet_bc_map.end() )
    {
      switch( bc_map_it->second )
	{
	  // No slip boundary condition
	case GRINS::NO_SLIP:
	  {
	    _bound_conds.apply_dirichlet( context, request_jacobian, _u_var, 0.0 );
	    
	    _bound_conds.apply_dirichlet( context, request_jacobian, _v_var, 0.0 );

	    if( _dim == 3 )
	      _bound_conds.apply_dirichlet( context, request_jacobian, _w_var, 0.0 );
	  }
	  break;

	  // Prescribed constant velocity
	case GRINS::PRESCRIBED_VELOCITY:
	  {
	    _bound_conds.apply_dirichlet( context, request_jacobian, 
					  _u_var, _vel_boundary_values[boundary_id][0] );

	    _bound_conds.apply_dirichlet( context, request_jacobian, 
					  _v_var, _vel_boundary_values[boundary_id][1] );
	    
	    if( _dim == 3 )
	      _bound_conds.apply_dirichlet( context, request_jacobian, 
					    _w_var, _vel_boundary_values[boundary_id][2] );
	  }
	  break;

	  // Inflow 
	case GRINS::INFLOW:
	  {
            std::map< GRINS::VariableIndex,GRINS::DirichletFuncObj* >& bc_map = _dirichlet_bound_funcs[boundary_id];
	    
	    std::map< GRINS::VariableIndex,GRINS::DirichletFuncObj* >::iterator u_it = bc_map.find( _u_var );
	    std::map< GRINS::VariableIndex,GRINS::DirichletFuncObj* >::iterator v_it = bc_map.find( _v_var );
	    std::map< GRINS::VariableIndex,GRINS::DirichletFuncObj* >::iterator w_it;

	    if( _dim == 3 ) w_it = bc_map.find( _w_var ); 
            
            if( u_it == bc_map.end() )
	      {
		_bound_conds.apply_dirichlet( context, request_jacobian, _u_var, 0.0 );
	      }
            else
	      {
		_bound_conds.apply_dirichlet( context, request_jacobian, _u_var, u_it->second );
	      }
	    if( v_it == bc_map.end() )
	      {
		_bound_conds.apply_dirichlet( context, request_jacobian, _v_var, 0.0 );
	      }
            else
	      {
		_bound_conds.apply_dirichlet( context, request_jacobian, _v_var, v_it->second );
	      }
            if( _dim == 3 )
	      {
		if( w_it == bc_map.end() )
		  {
		    _bound_conds.apply_dirichlet( context, request_jacobian, _w_var, 0.0 );
		  }
		else
		  {
		    _bound_conds.apply_dirichlet( context, request_jacobian, _w_var, w_it->second );
		  }
	      }
	  }
	  break;

	default:
	  {
	    std::cerr << "Error: Invalid Dirichlet BC type for IncompressibleNavierStokes."
		      << std::endl;
	    libmesh_error();
	  }

	} // End switch on bc type
    } // End if statement

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("IncompressibleNavierStokes::side_constraint");
#endif

  return request_jacobian;
}

bool GRINS::IncompressibleNavierStokes::mass_residual( bool request_jacobian,
						       libMesh::DiffContext& context,
						       libMesh::FEMSystem* system )
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // Element Jacobian * quadrature weights for interior integration
  // We assume the same for each flow variable
  const std::vector<Real> &JxW = 
    c.element_fe_var[_u_var]->get_JxW();

  // The shape functions at interior quadrature points.
  // We assume the same for each flow variable
  const std::vector<std::vector<Real> >& u_phi = 
    c.element_fe_var[_u_var]->get_phi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[_u_var].size();

  // for convenience
  if (_dim != 3)
    _w_var = _u_var;

  // The subvectors and submatrices we need to fill:
  DenseSubVector<Real> &F_u = *c.elem_subresiduals[_u_var];
  DenseSubVector<Real> &F_v = *c.elem_subresiduals[_v_var];
  DenseSubVector<Real> &F_w = *c.elem_subresiduals[_w_var];

  DenseSubMatrix<Real> &M_uu = *c.elem_subjacobians[_u_var][_u_var];
  DenseSubMatrix<Real> &M_vv = *c.elem_subjacobians[_v_var][_v_var];
  DenseSubMatrix<Real> &M_ww = *c.elem_subjacobians[_w_var][_w_var];

  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp = 0; qp != n_qpoints; ++qp)
    {
      // For the mass residual, we need to be a little careful.
      // The time integrator is handling the time-discretization
      // for us so we need to supply M(u_fixed)*u for the residual.
      // u_fixed will be given by the fixed_interior_* functions
      // while u will be given by the interior_* functions.
      Real u_dot = c.interior_value(_u_var, qp);
      Real v_dot = c.interior_value(_v_var, qp);

      Real w_dot = 0.0;

      if( _dim == 3 )
	Real w_dot = c.interior_value(_w_var, qp);
      
      for (unsigned int i = 0; i != n_u_dofs; ++i)
        {
	  F_u(i) += JxW[qp]*_rho*u_dot*u_phi[i][qp];
	  F_v(i) += JxW[qp]*_rho*v_dot*u_phi[i][qp];

	  if( _dim == 3 )
	    F_w(i) += JxW[qp]*_rho*w_dot*u_phi[i][qp];
	  
	  if( request_jacobian )
              {
		for (unsigned int j=0; j != n_u_dofs; j++)
		  {
		    // Assuming rho is constant w.r.t. u, v, w
		    // and T (if Boussinesq added).
		    Real value = JxW[qp]*_rho*u_phi[i][qp]*u_phi[j][qp];

		    M_uu(i,j) += value;
		    M_vv(i,j) += value;

		    if( _dim == 3)
		      {
			M_ww(i,j) += value;
		      }

		  } // End dof loop
	      } // End Jacobian check
	} // End dof loop
    } // End quadrature loop

  return request_jacobian;
}
