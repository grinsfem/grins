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
#include "grins/stokes.h"

// GRINS
#include "grins_config.h"

// libMesh
#include "libmesh/fem_context.h"
#include "libmesh/quadrature.h"

namespace GRINS
{

  Stokes::Stokes(const std::string& physics_name, const GetPot& input )
    : IncompressibleNavierStokesBase(physics_name,input),
      _p_pinning(input,physics_name),
      _pin_pressure( input("Physics/"+stokes+"/pin_pressure", false ) )
  {
    // This is deleted in the base class
    _bc_handler = new IncompressibleNavierStokesBCHandling( physics_name, input );

    return;
  }

  Stokes::~Stokes()
  {
    return;
  }

  void Stokes::element_time_derivative( bool compute_jacobian,
					libMesh::FEMContext& context,
					CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("Stokes::element_time_derivative");
#endif

    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.dof_indices_var[_u_var].size();
    const unsigned int n_p_dofs = context.dof_indices_var[_p_var].size();

    // Check number of dofs is same for _u_var, v_var and w_var.
    libmesh_assert (n_u_dofs == context.dof_indices_var[_v_var].size());
    if (_dim == 3)
      libmesh_assert (n_u_dofs == context.dof_indices_var[_w_var].size());

    // We get some references to cell-specific data that
    // will be used to assemble the linear system.

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.element_fe_var[_u_var]->get_JxW();

    // The velocity shape function gradients (in global coords.)
    // at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.element_fe_var[_u_var]->get_dphi();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& p_phi =
      context.element_fe_var[_p_var]->get_phi();

    // The subvectors and submatrices we need to fill:
    //
    // K_{\alpha \beta} = R_{\alpha},{\beta} = \partial{ R_{\alpha} } / \partial{ {\beta} } (where R denotes residual)
    // e.g., for \alpha = v and \beta = u we get: K{vu} = R_{v},{u}
    // Note that Kpu, Kpv, Kpw and Fp comes as constraint.
    //
    if (_dim != 3)
      _w_var = _u_var; // for convenience

    libMesh::DenseSubMatrix<libMesh::Number> &Kuu = *context.elem_subjacobians[_u_var][_u_var]; // R_{u},{u}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvv = *context.elem_subjacobians[_v_var][_v_var]; // R_{v},{v}
    libMesh::DenseSubMatrix<libMesh::Number> &Kww = *context.elem_subjacobians[_w_var][_w_var]; // R_{w},{w}

    libMesh::DenseSubMatrix<libMesh::Number> &Kup = *context.elem_subjacobians[_u_var][_p_var]; // R_{u},{p}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvp = *context.elem_subjacobians[_v_var][_p_var]; // R_{v},{p}
    libMesh::DenseSubMatrix<libMesh::Number> &Kwp = *context.elem_subjacobians[_w_var][_p_var]; // R_{w},{p}

    libMesh::DenseSubVector<libMesh::Number> &Fu = *context.elem_subresiduals[_u_var]; // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = *context.elem_subresiduals[_v_var]; // R_{v}
    libMesh::DenseSubVector<libMesh::Number> &Fw = *context.elem_subresiduals[_w_var]; // R_{w}

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
	libMesh::Number p, u, v, w;
	p = context.interior_value(_p_var, qp);
	u = context.interior_value(_u_var, qp);
	v = context.interior_value(_v_var, qp);
	if (_dim == 3)
	  w = context.interior_value(_w_var, qp);

	libMesh::Gradient grad_u, grad_v, grad_w;
	grad_u = context.interior_gradient(_u_var, qp);
	grad_v = context.interior_gradient(_v_var, qp);
	if (_dim == 3)
	  grad_w = context.interior_gradient(_w_var, qp);

	libMesh::NumberVectorValue Uvec (u,v);
	if (_dim == 3)
	  Uvec(2) = w;

	// First, an i-loop over the velocity degrees of freedom.
	// We know that n_u_dofs == n_v_dofs so we can compute contributions
	// for both at the same time.
	for (unsigned int i=0; i != n_u_dofs; i++)
	  {
	    Fu(i) += JxW[qp] *
	      ( p*u_gradphi[i][qp](0)              // pressure term
		-_mu*(u_gradphi[i][qp]*grad_u) ); // diffusion term

	    Fv(i) += JxW[qp] *
	      ( p*u_gradphi[i][qp](1)              // pressure term
		-_mu*(u_gradphi[i][qp]*grad_v) ); // diffusion term
	    if (_dim == 3)
	      {
		Fw(i) += JxW[qp] *
		  ( p*u_gradphi[i][qp](2)              // pressure term
		    -_mu*(u_gradphi[i][qp]*grad_w) ); // diffusion term
	      }

	    if (compute_jacobian)
	      {
		for (unsigned int j=0; j != n_u_dofs; j++)
		  {
		    // TODO: precompute some terms like:
		    //   (Uvec*vel_gblgradphivec[j][qp]),
		    //   vel_phi[i][qp]*vel_phi[j][qp],
		    //   (vel_gblgradphivec[i][qp]*vel_gblgradphivec[j][qp])

		    Kuu(i,j) += JxW[qp] *
		      (-_mu*(u_gradphi[i][qp]*u_gradphi[j][qp])); // diffusion term

		    Kvv(i,j) += JxW[qp] *
		      (-_mu*(u_gradphi[i][qp]*u_gradphi[j][qp])); // diffusion term

		    if (_dim == 3)
		      {
			Kww(i,j) += JxW[qp] *
			  (-_mu*(u_gradphi[i][qp]*u_gradphi[j][qp])); // diffusion term
		      }
		  } // end of the inner dof (j) loop

		// Matrix contributions for the up, vp and wp couplings
		for (unsigned int j=0; j != n_p_dofs; j++)
		  {
		    Kup(i,j) += JxW[qp]*u_gradphi[i][qp](0)*p_phi[j][qp];
		    Kvp(i,j) += JxW[qp]*u_gradphi[i][qp](1)*p_phi[j][qp];
		    if (_dim == 3)
		      Kwp(i,j) += JxW[qp]*u_gradphi[i][qp](2)*p_phi[j][qp];
		  } // end of the inner dof (j) loop

	      } // end - if (compute_jacobian && context.elem_solution_derivative)

	  } // end of the outer dof (i) loop
      } // end of the quadrature point (qp) loop

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("Stokes::element_time_derivative");
#endif

    return;
  }

  void Stokes::element_constraint( bool compute_jacobian,
				   libMesh::FEMContext& context,
				   CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("Stokes::element_constraint");
#endif

    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.dof_indices_var[_u_var].size();
    const unsigned int n_p_dofs = context.dof_indices_var[_p_var].size();

    // We get some references to cell-specific data that
    // will be used to assemble the linear system.

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.element_fe_var[_u_var]->get_JxW();

    // The velocity shape function gradients (in global coords.)
    // at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.element_fe_var[_u_var]->get_dphi();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& p_phi =
      context.element_fe_var[_p_var]->get_phi();

    // The subvectors and submatrices we need to fill:
    //
    // Kpu, Kpv, Kpw, Fp
    //
    if (_dim != 3)
      _w_var = _u_var; // for convenience

    libMesh::DenseSubMatrix<libMesh::Number> &Kpu = *context.elem_subjacobians[_p_var][_u_var]; // R_{p},{u}
    libMesh::DenseSubMatrix<libMesh::Number> &Kpv = *context.elem_subjacobians[_p_var][_v_var]; // R_{p},{v}
    libMesh::DenseSubMatrix<libMesh::Number> &Kpw = *context.elem_subjacobians[_p_var][_w_var]; // R_{p},{w}

    libMesh::DenseSubVector<libMesh::Number> &Fp = *context.elem_subresiduals[_p_var]; // R_{p}

    // Add the constraint given by the continuity equation.
    unsigned int n_qpoints = context.element_qrule->n_points();
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	// Compute the velocity gradient at the old Newton iterate.
	libMesh::Gradient grad_u, grad_v, grad_w;
	grad_u = context.interior_gradient(_u_var, qp);
	grad_v = context.interior_gradient(_v_var, qp);
	if (_dim == 3)
	  grad_w = context.interior_gradient(_w_var, qp);

	// Now a loop over the pressure degrees of freedom.  This
	// computes the contributions of the continuity equation.
	for (unsigned int i=0; i != n_p_dofs; i++)
	  {
	    Fp(i) += JxW[qp] * p_phi[i][qp] *
	      (grad_u(0) + grad_v(1));
	    if (_dim == 3)
	      Fp(i) += JxW[qp] * p_phi[i][qp] *
		(grad_w(2));

	    if (compute_jacobian)
	      {
		for (unsigned int j=0; j != n_u_dofs; j++)
		  {
		    Kpu(i,j) += JxW[qp]*p_phi[i][qp]*u_gradphi[j][qp](0);
		    Kpv(i,j) += JxW[qp]*p_phi[i][qp]*u_gradphi[j][qp](1);
		    if (_dim == 3)
		      Kpw(i,j) += JxW[qp]*p_phi[i][qp]*u_gradphi[j][qp](2);
		  } // end of the inner dof (j) loop

	      } // end - if (compute_jacobian && context.elem_solution_derivative)

	  } // end of the outer dof (i) loop
      } // end of the quadrature point (qp) loop


    // Pin p = p_value at p_point
    if( _pin_pressure )
      {
	_p_pinning.pin_value( context, compute_jacobian, _p_var );
      }
  

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("Stokes::element_constraint");
#endif

    return;
  }

  void Stokes::mass_residual( bool compute_jacobian,
			      libMesh::FEMContext& context,
			      CachedValues& /*cache*/)
  {
    // Element Jacobian * quadrature weights for interior integration
    // We assume the same for each flow variable
    const std::vector<libMesh::Real> &JxW = 
      context.element_fe_var[_u_var]->get_JxW();

    // The shape functions at interior quadrature points.
    // We assume the same for each flow variable
    const std::vector<std::vector<libMesh::Real> >& u_phi = 
      context.element_fe_var[_u_var]->get_phi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = context.dof_indices_var[_u_var].size();

    // for convenience
    if (_dim != 3)
      _w_var = _u_var;

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<libMesh::Real> &F_u = *context.elem_subresiduals[_u_var];
    libMesh::DenseSubVector<libMesh::Real> &F_v = *context.elem_subresiduals[_v_var];
    libMesh::DenseSubVector<libMesh::Real> &F_w = *context.elem_subresiduals[_w_var];

    libMesh::DenseSubMatrix<libMesh::Real> &M_uu = *context.elem_subjacobians[_u_var][_u_var];
    libMesh::DenseSubMatrix<libMesh::Real> &M_vv = *context.elem_subjacobians[_v_var][_v_var];
    libMesh::DenseSubMatrix<libMesh::Real> &M_ww = *context.elem_subjacobians[_w_var][_w_var];

    unsigned int n_qpoints = context.element_qrule->n_points();

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
	// For the mass residual, we need to be a little careful.
	// The time integrator is handling the time-discretization
	// for us so we need to supply M(u_fixed)*u for the residual.
	// u_fixed will be given by the fixed_interior_* functions
	// while u will be given by the interior_* functions.
	libMesh::Real u_dot = context.interior_value(_u_var, qp);
	libMesh::Real v_dot = context.interior_value(_v_var, qp);

	libMesh::Real w_dot = 0.0;

	if( _dim == 3 )
	  w_dot = context.interior_value(_w_var, qp);
      
	for (unsigned int i = 0; i != n_u_dofs; ++i)
	  {
	    F_u(i) += JxW[qp]*_rho*u_dot*u_phi[i][qp];
	    F_v(i) += JxW[qp]*_rho*v_dot*u_phi[i][qp];

	    if( _dim == 3 )
	      F_w(i) += JxW[qp]*_rho*w_dot*u_phi[i][qp];
	  
	    if( compute_jacobian )
              {
		for (unsigned int j=0; j != n_u_dofs; j++)
		  {
		    // Assuming rho is constant w.r.t. u, v, w
		    // and T (if Boussinesq added).
		    libMesh::Real value = JxW[qp]*_rho*u_phi[i][qp]*u_phi[j][qp];

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

    return;
  }

} // namespace GRINS
