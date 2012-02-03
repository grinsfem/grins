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

#include "axisym_heat_transfer.h"

void GRINS::AxisymmetricHeatTransfer::read_input_options( const GetPot& input )
{
  this->_T_FE_family =
    libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/AxisymmetricHeatTransfer/FE_family", "LAGRANGE") );

  this->_T_order =
    libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/AxisymmetricHeatTransfer/T_order", "SECOND") );

  this->_V_FE_family =
    libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/AxisymIncompNS/FE_family", "LAGRANGE") );

  this->_V_order =
    libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/AxisymIncompNS/V_order", "SECOND") );

  this->_rho = input("Physics/AxisymmetricHeatTransfer/rho", 1.0); //TODO: same as Incompressible NS
  this->_Cp  = input("Physics/AxisymmetricHeatTransfer/Cp", 1.0);
  this->_k  = input("Physics/AxisymmetricHeatTransfer/k", 1.0);

  this->_T_var_name = input("Physics/VariableNames/Temperature", GRINS::T_var_name_default );

  // registered/non-owned variable names
  this->_u_r_var_name = input("Physics/VariableNames/r_velocity", GRINS::u_r_var_name_default );
  this->_u_z_var_name = input("Physics/VariableNames/z_velocity", GRINS::u_z_var_name_default );
  

  // Read boundary condition info
  /** \todo We ought to be able to put this in the base class somehow so
            that it doesn't have to be rewritten for every physics class.
	    Then, the physics only handles the specifics, e.g. reading
	    in boundary velocities. */
  int num_ids = input.vector_variable_size("Physics/AxisymmetricHeatTransfer/bc_ids");
  int num_bcs = input.vector_variable_size("Physics/AxisymmetricHeatTransfer/bc_types");

  if( num_ids != num_bcs )
    {
      std::cerr << "Error: Must specify equal number of boundary ids and boundary conditions"
		<< std::endl;
      libmesh_error();
    }

  for( int i = 0; i < num_ids; i++ )
    {
      int bc_id = input("Physics/AxisymmetricHeatTransfer/bc_ids", -1, i );
      std::string bc_type_in = input("Physics/AxisymmetricHeatTransfer/bc_types", "NULL", i );

      GRINS::BC_TYPES bc_type = _bound_conds.string_to_enum( bc_type_in );

      std::stringstream ss;
      ss << bc_id;
      std::string bc_id_string = ss.str();

      // Now read in auxillary boundary condition information
      switch(bc_type)
	{
	case GRINS::ISOTHERMAL_WALL:
	  {
	    _dirichlet_bc_map[bc_id] = bc_type;

	    _T_boundary_values[bc_id] = 
	      input("Physics/AxisymmetricHeatTransfer/T_wall_"+bc_id_string, 0.0 );
	  }
	  break;

	case GRINS::ADIABATIC_WALL:
	  {
	    _neumann_bc_map[bc_id] = bc_type;
	  }
	  break;

	case GRINS::PRESCRIBED_HEAT_FLUX:
	  {
	    _neumann_bc_map[bc_id] = bc_type;

	    libMesh::Point q_in;
	    
	    int num_q_components = input.vector_variable_size("Physics/AxisymmetricHeatTransfer/q_wall_"+bc_id_string);
	    
	    if( num_q_components > 2 )
	      {
		std::cerr << "Error: Cannot have more than 2 components of heat flux for the axisymmetric case."
			  << std::endl
			  << "       Boundary id " << bc_id << " contains "
			  << num_q_components << "."
			  << std::endl;
		libmesh_error();
	      }

	    for( int i = 0; i < num_q_components; i++ )
	      {
		q_in(i) = input("Physics/AxisymmetricHeatTransfer/q_wall_"+bc_id_string, 0.0, i );
	      }
	      _q_boundary_values[bc_id] = q_in;
	  }
	  break;
	case GRINS::GENERAL_HEAT_FLUX:
	  {
	    _neumann_bc_map[bc_id] = bc_type;
	  }
	  break;
	case GRINS::AXISYMMETRIC:
	  {
	    _neumann_bc_map[bc_id] = bc_type;
	  }
	  break;

	default:
	  {
	    std::cerr << "Error: Invalid boundary condition type for AxisymmetricHeatTransfer."
		      << std::endl;
	    libmesh_error();
	  }

	}// End switch(bc_type)
    } // End loop on bc_id

  return;
}

void GRINS::AxisymmetricHeatTransfer::init_variables( libMesh::FEMSystem* system )
{
  // Get libMesh to assign an index for each variable
  this->_dim = system->get_mesh().mesh_dimension();

  _T_var = system->add_variable( _T_var_name, _T_order, _T_FE_family);

  // If these are already added, then we just get the index. 
  _u_r_var = system->add_variable(_u_r_var_name, _V_order, _V_FE_family);
  _u_z_var = system->add_variable(_u_z_var_name, _V_order, _V_FE_family);

  return;
}

void GRINS::AxisymmetricHeatTransfer::set_time_evolving_vars( libMesh::FEMSystem* system )
{
  // Tell the system to march temperature forward in time
  system->time_evolving(_T_var);
  return;
}

void GRINS::AxisymmetricHeatTransfer::init_context( libMesh::DiffContext &context )
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

  // _u_var is registered so can we assume things related to _u_var
  // are available in FEMContext

  return;
}

bool GRINS::AxisymmetricHeatTransfer::element_time_derivative( bool request_jacobian,
							       libMesh::DiffContext& context,
							       libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("AxisymmetricHeatTransfer::element_time_derivative");
#endif

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // The number of local degrees of freedom in each variable.
  const unsigned int n_T_dofs = c.dof_indices_var[_T_var].size();
  const unsigned int n_u_dofs = c.dof_indices_var[_u_r_var].size();

  //TODO: check n_T_dofs is same as n_u_dofs, n_v_dofs, n_w_dofs

  // We get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[_T_var]->get_JxW();

  // The temperature shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& T_phi =
    c.element_fe_var[_T_var]->get_phi();

  // The velocity shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& vel_phi =
    c.element_fe_var[_u_r_var]->get_phi();

  // The temperature shape function gradients (in global coords.)
  // at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
    c.element_fe_var[_T_var]->get_dphi();

  // Physical location of the quadrature points
  const std::vector<libMesh::Point>& u_qpoint =
    c.element_fe_var[_u_r_var]->get_xyz();

  // The subvectors and submatrices we need to fill:
  libMesh::DenseSubVector<Number> &FT = *c.elem_subresiduals[_T_var]; // R_{T}

  libMesh::DenseSubMatrix<Number> &KTT = *c.elem_subjacobians[_T_var][_T_var]; // R_{T},{T}

  libMesh::DenseSubMatrix<Number> &KTr = *c.elem_subjacobians[_T_var][_u_r_var]; // R_{T},{r}
  libMesh::DenseSubMatrix<Number> &KTz = *c.elem_subjacobians[_T_var][_u_z_var]; // R_{T},{z}


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
      libMesh::Number T, u_r, u_z;
      T = c.interior_value(_T_var, qp);
      u_r = c.interior_value(_u_r_var, qp);
      u_z = c.interior_value(_u_z_var, qp);

      libMesh::Gradient grad_T;
      grad_T = c.interior_gradient(_T_var, qp);

      libMesh::NumberVectorValue U (u_r,u_z);

      // First, an i-loop over the  degrees of freedom.
      for (unsigned int i=0; i != n_T_dofs; i++)
        {
          FT(i) += JxW[qp]*r*
                   (-_rho*_Cp*T_phi[i][qp]*(U*grad_T)    // convection term
                    -_k*(T_gradphi[i][qp]*grad_T) );  // diffusion term

          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);

              for (unsigned int j=0; j != n_T_dofs; j++)
                {
                  // TODO: precompute some terms like:
                  //   _rho*_Cp*T_phi[i][qp]*(vel_phi[j][qp]*T_grad_phi[j][qp])

                  KTT(i,j) += JxW[qp]*r*
                              (-_rho*_Cp*T_phi[i][qp]*(U*T_gradphi[j][qp])  // convection term
                               -_k*(T_gradphi[i][qp]*T_gradphi[j][qp])); // diffusion term
                } // end of the inner dof (j) loop

              // Matrix contributions for the Tu, Tv and Tw couplings (n_T_dofs same as n_u_dofs, n_v_dofs and n_w_dofs)
              for (unsigned int j=0; j != n_u_dofs; j++)
                {
                  KTr(i,j) += JxW[qp]*r*(-_rho*_Cp*T_phi[i][qp]*(vel_phi[j][qp]*grad_T(0)));
                  KTz(i,j) += JxW[qp]*r*(-_rho*_Cp*T_phi[i][qp]*(vel_phi[j][qp]*grad_T(1)));
                } // end of the inner dof (j) loop

            } // end - if (request_jacobian && c.elem_solution_derivative)

        } // end of the outer dof (i) loop
    } // end of the quadrature point (qp) loop

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("AxisymmetricHeatTransfer::element_time_derivative");
#endif

  return request_jacobian;
}

bool GRINS::AxisymmetricHeatTransfer::element_constraint( bool request_jacobian,
							    libMesh::DiffContext& context,
							    libMesh::FEMSystem* system )
{
  return request_jacobian;
}


bool GRINS::AxisymmetricHeatTransfer::side_time_derivative( bool request_jacobian,
							      libMesh::DiffContext& context,
							      libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("AxisymmetricHeatTransfer::side_time_derivative");
#endif

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  const GRINS::BoundaryID boundary_id =
    system->get_mesh().boundary_info->boundary_id(c.elem, c.side);
  libmesh_assert (boundary_id != libMesh::BoundaryInfo::invalid_id);

  std::map< GRINS::BoundaryID, GRINS::BC_TYPES>::const_iterator 
    bc_map_it = _neumann_bc_map.find( boundary_id );

   /* We assume that if you didn't put a boundary id in, then you didn't want to
     set a boundary condition on that boundary. */
  if( bc_map_it != _neumann_bc_map.end() )
    {
      switch( bc_map_it->second )
	{
	  // Zero heat flux
	case GRINS::ADIABATIC_WALL:
	  // Don't need to do anything: q = 0 in this case
	  break;

	  // Prescribed constant heat flux
	case GRINS::PRESCRIBED_HEAT_FLUX:
	  {
	    _bound_conds.apply_neumann( context, _T_var, -1.0,
					_q_boundary_values[boundary_id] );
	  }
	  break;
	case GRINS::GENERAL_HEAT_FLUX:
	  {
	    GRINS::NeumannBCsMap& bc_map = _neumann_bound_funcs[boundary_id];
	    
	    GRINS::NeumannBCsMap::iterator T_it = bc_map.find( _T_var );
	    
	    _bound_conds.apply_neumann( context, request_jacobian, _T_var, -1.0, T_it->second );
	  }
	  break;
	case GRINS::AXISYMMETRIC:
	  // Don't need to do anything for dT/dr = 0
	  break;
	default:
	  {
	    std::cerr << "Error: Invalid Neumann BC type for AxisymmetricHeatTransfer."
		      << std::endl;
	    libmesh_error();
	  }
	} // End switch
    } // End if

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("AxisymmetricHeatTransfer::side_time_derivative");
#endif

  return request_jacobian;
}

bool GRINS::AxisymmetricHeatTransfer::side_constraint( bool request_jacobian,
							 libMesh::DiffContext& context,
							 libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("AxisymmetricHeatTransfer::side_constraint");
#endif

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

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
	  // Prescribed constant temperature
	case GRINS::ISOTHERMAL_WALL:
	  {
	    _bound_conds.apply_dirichlet( context, request_jacobian,
					  _T_var, _T_boundary_values[boundary_id] );
	  }
	  break;
	default:
	  {
	    std::cerr << "Error: Invalid Dirichlet BC type for AxisymmetricHeatTransfer."
		      << std::endl;
	    libmesh_error();
	  }

	} // End switch
    } // End if

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("AxisymmetricHeatTransfer::side_constraint");
#endif

  return request_jacobian;
}

bool GRINS::AxisymmetricHeatTransfer::mass_residual( bool request_jacobian,
						     libMesh::DiffContext& context,
						     libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("AxisymmetricHeatTransfer::mass_residual");
#endif

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = 
    c.element_fe_var[_T_var]->get_JxW();

  // The shape functions at interior quadrature points.
  const std::vector<std::vector<Real> >& phi = 
    c.element_fe_var[_T_var]->get_phi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_T_dofs = c.dof_indices_var[_T_var].size();

  // Physical location of the quadrature points
  const std::vector<libMesh::Point>& u_qpoint =
    c.element_fe_var[_u_r_var]->get_xyz();

  // The subvectors and submatrices we need to fill:
  DenseSubVector<Real> &F = *c.elem_subresiduals[_T_var];

  DenseSubMatrix<Real> &M = *c.elem_subjacobians[_T_var][_T_var];

  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp = 0; qp != n_qpoints; ++qp)
    {
      const libMesh::Number r = u_qpoint[qp](0);

      // For the mass residual, we need to be a little careful.
      // The time integrator is handling the time-discretization
      // for us so we need to supply M(u_fixed)*u for the residual.
      // u_fixed will be given by the fixed_interior_* functions
      // while u will be given by the interior_* functions.
      Real T_dot = c.interior_value(_T_var, qp);

      for (unsigned int i = 0; i != n_T_dofs; ++i)
        {
          F(i) += JxW[qp]*r*(_rho*_Cp*T_dot*phi[i][qp] );

          if( request_jacobian )
              {
                for (unsigned int j=0; j != n_T_dofs; j++)
                  {
		    // We're assuming rho, cp are constant w.r.t. T here.
                    M(i,j) += JxW[qp]*r*_rho*_Cp*phi[j][qp]*phi[i][qp] ;
                  }
              }// End of check on Jacobian
          
        } // End of element dof loop
      
    } // End of the quadrature point loop

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("AxisymmetricHeatTransfer::mass_residual");
#endif

  return request_jacobian;
}
