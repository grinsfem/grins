//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "grins/axisym_heat_transfer.h"

// GRINS
#include "grins/axisym_heat_transfer_bc_handling.h"
#include "grins/constant_conductivity.h"

// libMesh
#include "libmesh/utility.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/fem_context.h"
#include "libmesh/quadrature.h"

namespace GRINS
{

  template< class Conductivity>
  AxisymmetricHeatTransfer<Conductivity>::AxisymmetricHeatTransfer( const std::string& physics_name,
								    const GetPot& input)
    : Physics(physics_name, input)
  {
    this->read_input_options(input);
  
    // This is deleted in the base class
    _bc_handler = new AxisymmetricHeatTransferBCHandling( physics_name, input );

    return;
  }

  template< class Conductivity>
  AxisymmetricHeatTransfer<Conductivity>::~AxisymmetricHeatTransfer( )
  {
    return;
  }

  template< class Conductivity>
  void AxisymmetricHeatTransfer<Conductivity>::read_input_options( const GetPot& input )
  {
    this->_T_FE_family =
      libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+axisymmetric_heat_transfer+"/FE_family", "LAGRANGE") );

    this->_T_order =
      libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+axisymmetric_heat_transfer+"/T_order", "SECOND") );

    this->_V_FE_family =
      libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+axisymmetric_incomp_navier_stokes+"/FE_family", "LAGRANGE") );

    this->_V_order =
      libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+axisymmetric_incomp_navier_stokes+"/V_order", "SECOND") );

    this->_rho = input("Physics/"+axisymmetric_heat_transfer+"/rho", 1.0); //TODO: same as Incompressible NS
    this->_Cp  = input("Physics/"+axisymmetric_heat_transfer+"/Cp", 1.0);

    this->_k.read_input_options( input );

    this->_T_var_name = input("Physics/VariableNames/Temperature", T_var_name_default );

    // registered/non-owned variable names
    this->_u_r_var_name = input("Physics/VariableNames/r_velocity", u_r_var_name_default );
    this->_u_z_var_name = input("Physics/VariableNames/z_velocity", u_z_var_name_default );

    return;
  }

  template< class Conductivity>
  void AxisymmetricHeatTransfer<Conductivity>::init_variables( libMesh::FEMSystem* system )
  {
    // Get libMesh to assign an index for each variable
    this->_dim = system->get_mesh().mesh_dimension();

    _T_var = system->add_variable( _T_var_name, _T_order, _T_FE_family);

    // If these are already added, then we just get the index. 
    _u_r_var = system->add_variable(_u_r_var_name, _V_order, _V_FE_family);
    _u_z_var = system->add_variable(_u_z_var_name, _V_order, _V_FE_family);

    return;
  }

  template< class Conductivity>
  void AxisymmetricHeatTransfer<Conductivity>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Tell the system to march temperature forward in time
    system->time_evolving(_T_var);
    return;
  }

  template< class Conductivity>
  void AxisymmetricHeatTransfer<Conductivity>::init_context( libMesh::FEMContext& context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    context.element_fe_var[_T_var]->get_JxW();
    context.element_fe_var[_T_var]->get_phi();
    context.element_fe_var[_T_var]->get_dphi();
    context.element_fe_var[_T_var]->get_xyz();

    context.side_fe_var[_T_var]->get_JxW();
    context.side_fe_var[_T_var]->get_phi();
    context.side_fe_var[_T_var]->get_dphi();
    context.side_fe_var[_T_var]->get_xyz();

    // _u_var is registered so can we assume things related to _u_var
    // are available in FEMContext

    return;
  }

  template< class Conductivity>
  void AxisymmetricHeatTransfer<Conductivity>::element_time_derivative( bool compute_jacobian,
									libMesh::FEMContext& context,
									CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("AxisymmetricHeatTransfer::element_time_derivative");
#endif

    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.dof_indices_var[_T_var].size();
    const unsigned int n_u_dofs = context.dof_indices_var[_u_r_var].size();

    //TODO: check n_T_dofs is same as n_u_dofs, n_v_dofs, n_w_dofs

    // We get some references to cell-specific data that
    // will be used to assemble the linear system.

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.element_fe_var[_T_var]->get_JxW();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.element_fe_var[_T_var]->get_phi();

    // The velocity shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& vel_phi =
      context.element_fe_var[_u_r_var]->get_phi();

    // The temperature shape function gradients (in global coords.)
    // at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.element_fe_var[_T_var]->get_dphi();

    // Physical location of the quadrature points
    const std::vector<libMesh::Point>& u_qpoint =
      context.element_fe_var[_u_r_var]->get_xyz();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<Number> &FT = *context.elem_subresiduals[_T_var]; // R_{T}

    libMesh::DenseSubMatrix<Number> &KTT = *context.elem_subjacobians[_T_var][_T_var]; // R_{T},{T}

    libMesh::DenseSubMatrix<Number> &KTr = *context.elem_subjacobians[_T_var][_u_r_var]; // R_{T},{r}
    libMesh::DenseSubMatrix<Number> &KTz = *context.elem_subjacobians[_T_var][_u_z_var]; // R_{T},{z}


    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.element_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	const libMesh::Number r = u_qpoint[qp](0);
      
	// Compute the solution & its gradient at the old Newton iterate.
	libMesh::Number T, u_r, u_z;
	T = context.interior_value(_T_var, qp);
	u_r = context.interior_value(_u_r_var, qp);
	u_z = context.interior_value(_u_z_var, qp);

	libMesh::Gradient grad_T;
	grad_T = context.interior_gradient(_T_var, qp);

	libMesh::NumberVectorValue U (u_r,u_z);

	libMesh::Number k = this->_k( T );
	libMesh::Number dk_dT = this->_k.deriv( T );

	// First, an i-loop over the  degrees of freedom.
	for (unsigned int i=0; i != n_T_dofs; i++)
	  {
	    FT(i) += JxW[qp]*r*
	      (-_rho*_Cp*T_phi[i][qp]*(U*grad_T)    // convection term
	       -k*(T_gradphi[i][qp]*grad_T) );  // diffusion term

	    if (compute_jacobian && context.elem_solution_derivative)
	      {
		libmesh_assert (context.elem_solution_derivative == 1.0);

		for (unsigned int j=0; j != n_T_dofs; j++)
		  {
		    // TODO: precompute some terms like:
		    //   _rho*_Cp*T_phi[i][qp]*(vel_phi[j][qp]*T_grad_phi[j][qp])

		    KTT(i,j) += JxW[qp]*r*
		      (-_rho*_Cp*T_phi[i][qp]*(U*T_gradphi[j][qp])  // convection term
		       -k*(T_gradphi[i][qp]*T_gradphi[j][qp])); // diffusion term
		  } // end of the inner dof (j) loop

		//if( dk_dT != 0.0 )
		{
		  for (unsigned int j=0; j != n_T_dofs; j++)
		    {
		      // TODO: precompute some terms like:
		      KTT(i,j) -= JxW[qp]*r*( dk_dT*T_phi[j][qp]*T_gradphi[i][qp]*grad_T );
		    }
		}

		// Matrix contributions for the Tu, Tv and Tw couplings (n_T_dofs same as n_u_dofs, n_v_dofs and n_w_dofs)
		for (unsigned int j=0; j != n_u_dofs; j++)
		  {
		    KTr(i,j) += JxW[qp]*r*(-_rho*_Cp*T_phi[i][qp]*(vel_phi[j][qp]*grad_T(0)));
		    KTz(i,j) += JxW[qp]*r*(-_rho*_Cp*T_phi[i][qp]*(vel_phi[j][qp]*grad_T(1)));
		  } // end of the inner dof (j) loop

	      } // end - if (compute_jacobian && context.elem_solution_derivative)

	  } // end of the outer dof (i) loop
      } // end of the quadrature point (qp) loop

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("AxisymmetricHeatTransfer::element_time_derivative");
#endif

    return;
  }

  template< class Conductivity>
  void AxisymmetricHeatTransfer<Conductivity>::side_time_derivative( bool compute_jacobian,
								     libMesh::FEMContext& context,
								     CachedValues& cache )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("AxisymmetricHeatTransfer::side_time_derivative");
#endif

    std::vector<BoundaryID> ids = context.side_boundary_ids();

    for( std::vector<BoundaryID>::const_iterator it = ids.begin();
	 it != ids.end(); it++ )
      {
	libmesh_assert (*it != libMesh::BoundaryInfo::invalid_id);

	_bc_handler->apply_neumann_bcs( context, cache, compute_jacobian, *it );
      }

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("AxisymmetricHeatTransfer::side_time_derivative");
#endif

    return;
  }

  template< class Conductivity>
  void AxisymmetricHeatTransfer<Conductivity>::mass_residual( bool compute_jacobian,
							      libMesh::FEMContext& context,
							      CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("AxisymmetricHeatTransfer::mass_residual");
#endif

    // First we get some references to cell-specific data that
    // will be used to assemble the linear system.

    // Element Jacobian * quadrature weights for interior integration
    const std::vector<Real> &JxW = 
      context.element_fe_var[_T_var]->get_JxW();

    // The shape functions at interior quadrature points.
    const std::vector<std::vector<Real> >& phi = 
      context.element_fe_var[_T_var]->get_phi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_T_dofs = context.dof_indices_var[_T_var].size();

    // Physical location of the quadrature points
    const std::vector<libMesh::Point>& u_qpoint =
      context.element_fe_var[_u_r_var]->get_xyz();

    // The subvectors and submatrices we need to fill:
    DenseSubVector<Real> &F = *context.elem_subresiduals[_T_var];

    DenseSubMatrix<Real> &M = *context.elem_subjacobians[_T_var][_T_var];

    unsigned int n_qpoints = context.element_qrule->n_points();

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
	const libMesh::Number r = u_qpoint[qp](0);

	// For the mass residual, we need to be a little careful.
	// The time integrator is handling the time-discretization
	// for us so we need to supply M(u_fixed)*u for the residual.
	// u_fixed will be given by the fixed_interior_* functions
	// while u will be given by the interior_* functions.
	Real T_dot = context.interior_value(_T_var, qp);

	for (unsigned int i = 0; i != n_T_dofs; ++i)
	  {
	    F(i) += JxW[qp]*r*(_rho*_Cp*T_dot*phi[i][qp] );

	    if( compute_jacobian )
              {
                for (unsigned int j=0; j != n_T_dofs; j++)
                  {
		    // We're assuming rho, cp are constant w.r.t. T here.
                    M(i,j) += JxW[qp]*r*_rho*_Cp*phi[j][qp]*phi[i][qp] ;
                  }
              }// End of check on Jacobian
          
	  } // End of element dof loop
      
      } // End of the quadrature point loop

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("AxisymmetricHeatTransfer::mass_residual");
#endif

    return;
  }

} // namespace GRINS

template class GRINS::AxisymmetricHeatTransfer<GRINS::ConstantConductivity>;
