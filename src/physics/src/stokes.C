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

#include "stokes.h"

GRINS::Stokes::Stokes(const std::string& physics_name, const GetPot& input )
  : GRINS::IncompressibleNavierStokesBase(physics_name,input),
    _p_pinning(input,physics_name),
    _pin_pressure( input("Physics/"+stokes+"/pin_pressure", false ) )
{
  this->read_input_options(input);

  // This is deleted in the base class
  _bc_handler = new GRINS::IncompressibleNavierStokesBCHandling( physics_name, input );

  return;
}

GRINS::Stokes::~Stokes()
{
  return;
}

void GRINS::Stokes::read_input_options( const GetPot& input )
{
  return;
}

bool GRINS::Stokes::element_time_derivative( bool request_jacobian,
								 libMesh::DiffContext& context,
								 libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("Stokes::element_time_derivative");
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
                   ( p*vel_gblgradphivec[i][qp](0)              // pressure term
                    -_mu*(vel_gblgradphivec[i][qp]*graduvec) ); // diffusion term

          Fv(i) += JxW[qp] *
                   ( p*vel_gblgradphivec[i][qp](1)              // pressure term
                    -_mu*(vel_gblgradphivec[i][qp]*gradvvec) ); // diffusion term
          if (_dim == 3)
            {
              Fw(i) += JxW[qp] *
                       ( p*vel_gblgradphivec[i][qp](2)              // pressure term
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
                              (-_mu*(vel_gblgradphivec[i][qp]*vel_gblgradphivec[j][qp])); // diffusion term

                  Kvv(i,j) += JxW[qp] *
                              (-_mu*(vel_gblgradphivec[i][qp]*vel_gblgradphivec[j][qp])); // diffusion term

                  if (_dim == 3)
                    {
                      Kww(i,j) += JxW[qp] *
                                  (-_mu*(vel_gblgradphivec[i][qp]*vel_gblgradphivec[j][qp])); // diffusion term
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
  this->_timer->EndTimer("Stokes::element_time_derivative");
#endif

  return request_jacobian;
}

bool GRINS::Stokes::element_constraint( bool request_jacobian,
							    libMesh::DiffContext& context,
							    libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("Stokes::element_constraint");
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
      _p_pinning.pin_value( context, request_jacobian, _p_var );
    }
  

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("Stokes::element_constraint");
#endif

  return request_jacobian;
}


bool GRINS::Stokes::side_time_derivative( bool request_jacobian,
							      libMesh::DiffContext& context,
							      libMesh::FEMSystem* system )
{
  return request_jacobian;
}

bool GRINS::Stokes::side_constraint( bool request_jacobian,
							 libMesh::DiffContext& context,
							 libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  //this->_timer->BeginTimer("Stokes::side_constraint");
#endif

  //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

#ifdef USE_GRVY_TIMERS
  //this->_timer->EndTimer("Stokes::side_constraint");
#endif

  return request_jacobian;
}

bool GRINS::Stokes::mass_residual( bool request_jacobian,
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
