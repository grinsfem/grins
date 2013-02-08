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
#include "grins/heat_conduction.h"

// GRINS
#include "grins/heat_transfer_bc_handling.h"

// libMesh
#include "libmesh/string_to_enum.h"
#include "libmesh/quadrature.h"
#include "libmesh/fem_system.h"
#include "libmesh/fem_context.h"

namespace GRINS
{

  HeatConduction::HeatConduction( const GRINS::PhysicsName& physics_name, const GetPot& input )
    : Physics(physics_name,input)
  {
     this->_T_FE_family =
      libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/HeatConduction/FE_family", "LAGRANGE") );

    this->_T_order =
      libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/HeatConduction/T_order", "FIRST") );

    this->_rho = input("Physics/HeatConduction/rho", 1.0); //TODO: same as Incompressible NS
    this->_Cp  = input("Physics/HeatConduction/Cp", 1.0);
    this->_k  = input("Physics/HeatConduction/k", 1.0);

    this->_T_var_name = input("Physics/VariableNames/Temperature", T_var_name_default );

    // This is deleted in the base class
    _bc_handler = new HeatTransferBCHandling( physics_name, input );

    return;
  }

  HeatConduction::~HeatConduction()
  {
    return;
  }

  void HeatConduction::init_variables( libMesh::FEMSystem* system )
  {
    // Get libMesh to assign an index for each variable
    this->_dim = system->get_mesh().mesh_dimension();

    _T_var = system->add_variable( _T_var_name, this->_T_order, _T_FE_family);
    
    return;
  }

  void HeatConduction::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Tell the system to march temperature forward in time
    system->time_evolving(_T_var);

    return;
  }

  void HeatConduction::init_context( libMesh::FEMContext& context )
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

    return;
  }

  void HeatConduction::element_time_derivative( bool compute_jacobian,
						libMesh::FEMContext& context,
						CachedValues& /*cache*/ )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.dof_indices_var[_T_var].size();

    // We get some references to cell-specific data that
    // will be used to assemble the linear system.

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.element_fe_var[_T_var]->get_JxW();

    // The temperature shape function gradients (in global coords.)
    // at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.element_fe_var[_T_var]->get_dphi();

    const std::vector<libMesh::Point>& q_points = 
      context.element_fe_var[_T_var]->get_xyz();

    // The subvectors and submatrices we need to fill:
    //
    // K_{\alpha \beta} = R_{\alpha},{\beta} = \partial{ R_{\alpha} } / \partial{ {\beta} } (where R denotes residual)
    // e.g., for \alpha = T and \beta = v we get: K_{Tu} = R_{T},{u}
    //

    libMesh::DenseSubMatrix<Number> &KTT = *context.elem_subjacobians[_T_var][_T_var]; // R_{T},{T}

    libMesh::DenseSubVector<Number> &FT = *context.elem_subresiduals[_T_var]; // R_{T}

    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.element_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	// Compute the solution & its gradient at the old Newton iterate.
	libMesh::Gradient grad_T;
	grad_T = context.interior_gradient(_T_var, qp);

	const Real f = this->forcing( q_points[qp] );
	
	// First, an i-loop over the  degrees of freedom.
	for (unsigned int i=0; i != n_T_dofs; i++)
	  {
	    FT(i) += JxW[qp] *
	      (-_k*(T_gradphi[i][qp]*grad_T) + f);  // diffusion term

	    if (compute_jacobian)
	      {
		for (unsigned int j=0; j != n_T_dofs; j++)
		  {
		    // TODO: precompute some terms like:
		    //   _rho*_Cp*T_phi[i][qp]*(vel_phi[j][qp]*T_grad_phi[j][qp])

		    KTT(i,j) += JxW[qp] *
		      ( -_k*(T_gradphi[i][qp]*T_gradphi[j][qp]) ); // diffusion term
		  } // end of the inner dof (j) loop

	      } // end - if (compute_jacobian && context.elem_solution_derivative)

	  } // end of the outer dof (i) loop
      } // end of the quadrature point (qp) loop

    return;
  }

  void HeatConduction::mass_residual( bool compute_jacobian,
				      libMesh::FEMContext& context,
				      CachedValues& /*cache*/ )
  {
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

    // The subvectors and submatrices we need to fill:
    DenseSubVector<Real> &F = *context.elem_subresiduals[_T_var];

    DenseSubMatrix<Real> &M = *context.elem_subjacobians[_T_var][_T_var];

    unsigned int n_qpoints = context.element_qrule->n_points();

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
	// For the mass residual, we need to be a little careful.
	// The time integrator is handling the time-discretization
	// for us so we need to supply M(u_fixed)*u for the residual.
	// u_fixed will be given by the fixed_interior_* functions
	// while u will be given by the interior_* functions.
	Real T_dot = context.interior_value(_T_var, qp);

	for (unsigned int i = 0; i != n_T_dofs; ++i)
	  {
	    F(i) += JxW[qp]*(_rho*_Cp*T_dot*phi[i][qp] );

	    if( compute_jacobian )
              {
                for (unsigned int j=0; j != n_T_dofs; j++)
                  {
		    // We're assuming rho, cp are constant w.r.t. T here.
                    M(i,j) += JxW[qp]*_rho*_Cp*phi[j][qp]*phi[i][qp] ;
                  }
              }// End of check on Jacobian
          
	  } // End of element dof loop
      
      } // End of the quadrature point loop

    return;
  }

  inline
  Real HeatConduction::forcing( const libMesh::Point& p )
  {
    const Real x = p(0);
    const Real y = p(1);
    const Real z = p(2);

    return std::cos(.5*pi*x)*std::sin(.5*pi*y)*std::cos(.5*pi*z);
  }

} // namespace GRINS
