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

#include "low_mach_navier_stokes.h"

template<class Mu, class SH, class TC>
GRINS::LowMachNavierStokes<Mu,SH,TC>::LowMachNavierStokes(const std::string& physics_name, const GetPot& input)
  : GRINS::LowMachNavierStokesBase<Mu,SH,TC>(physics_name,input),
    _p_pinning(input,physics_name)
{
  this->read_input_options(input);

  // This is deleted in the base class
  this->_bc_handler = new GRINS::LowMachNavierStokesBCHandling( physics_name, input );

  return;
}

template<class Mu, class SH, class TC>
GRINS::LowMachNavierStokes<Mu,SH,TC>::~LowMachNavierStokes()
{
  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokes<Mu,SH,TC>::read_input_options( const GetPot& input )
{
  // Other quantities read in base class

  // Read pressure pinning information
  this->_pin_pressure = input("Physics/"+low_mach_navier_stokes+"/pin_pressure", false );
  
  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokes<Mu,SH,TC>::init_context( libMesh::DiffContext &context )
{
  // First call base class
  GRINS::LowMachNavierStokesBase<Mu,SH,TC>::init_context(context);

  libMesh::FEMContext &c = libmesh_cast_ref<libMesh::FEMContext&>(context);

  // We also need the side shape functions, etc.
  c.side_fe_var[this->_u_var]->get_JxW();
  c.side_fe_var[this->_u_var]->get_phi();
  c.side_fe_var[this->_u_var]->get_dphi();
  c.side_fe_var[this->_u_var]->get_xyz();

  c.side_fe_var[this->_T_var]->get_JxW();
  c.side_fe_var[this->_T_var]->get_phi();
  c.side_fe_var[this->_T_var]->get_dphi();
  c.side_fe_var[this->_T_var]->get_xyz();

  return;
}


template<class Mu, class SH, class TC>
bool GRINS::LowMachNavierStokes<Mu,SH,TC>::element_time_derivative( bool request_jacobian,
								    libMesh::DiffContext& context,
								    libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("LowMachNavierStokes::element_time_derivative");
#endif

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  this->assemble_mass_time_deriv( request_jacobian, c, system );
  this->assemble_momentum_time_deriv( request_jacobian, c, system );
  this->assemble_energy_time_deriv( request_jacobian, c, system );

  // Pin p = p_value at p_point
  if( this->_pin_pressure )
    {
      this->_p_pinning.pin_value( context, request_jacobian, this->_p_var);
    }

  if( this->_enable_thermo_press_calc )
    this->assemble_thermo_press_elem_time_deriv( request_jacobian, c, system );

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("LowMachNavierStokes::element_time_derivative");
#endif

  return request_jacobian;
}

template<class Mu, class SH, class TC>
bool GRINS::LowMachNavierStokes<Mu,SH,TC>::element_constraint( bool request_jacobian,
							       libMesh::DiffContext& context,
							       libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  //this->_timer->BeginTimer("LowMachNavierStokes::element_constraint");
#endif

  //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

#ifdef USE_GRVY_TIMERS
  //this->_timer->EndTimer("LowMachNavierStokes::element_constraint");
#endif

  return request_jacobian;
}


template<class Mu, class SH, class TC>
bool GRINS::LowMachNavierStokes<Mu,SH,TC>::side_time_derivative( bool request_jacobian,
								 libMesh::DiffContext& context,
								 libMesh::FEMSystem* system )
{
  if( this->_enable_thermo_press_calc )
    {
#ifdef USE_GRVY_TIMERS
      this->_timer->BeginTimer("LowMachNavierStokes::side_time_derivative");
#endif
      FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

      this->assemble_thermo_press_side_time_deriv( request_jacobian, c, system );

#ifdef USE_GRVY_TIMERS
      this->_timer->EndTimer("LowMachNavierStokes::side_time_derivative");
#endif
    }

  return request_jacobian;
}

template<class Mu, class SH, class TC>
bool GRINS::LowMachNavierStokes<Mu,SH,TC>::side_constraint( bool request_jacobian,
							    libMesh::DiffContext& context,
							    libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  //this->_timer->BeginTimer("LowMachNavierStokes::side_constraint");
#endif

  //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

#ifdef USE_GRVY_TIMERS
  //this->_timer->EndTimer("LowMachNavierStokes::side_constraint");
#endif

  return request_jacobian;
}

template<class Mu, class SH, class TC>
bool GRINS::LowMachNavierStokes<Mu,SH,TC>::mass_residual( bool request_jacobian,
							  libMesh::DiffContext& context,
							  libMesh::FEMSystem* system )
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  this->assemble_continuity_mass_residual( request_jacobian, c, system );

  this->assemble_momentum_mass_residual( request_jacobian, c, system );

  this->assemble_energy_mass_residual( request_jacobian, c, system );

  if( this->_enable_thermo_press_calc )
    this->assemble_thermo_press_mass_residual( request_jacobian, c, system );

  return request_jacobian;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokes<Mu,SH,TC>::assemble_mass_time_deriv( bool request_jacobian, 
								     libMesh::FEMContext& c, 
								     libMesh::FEMSystem* system )
{
  // The number of local degrees of freedom in each variable.
  const unsigned int n_p_dofs = c.dof_indices_var[this->_p_var].size();

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[this->_u_var]->get_JxW();

  // The pressure shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& p_phi =
    c.element_fe_var[this->_p_var]->get_phi();

  libMesh::DenseSubVector<Number> &Fp = *c.elem_subresiduals[this->_p_var]; // R_{p}

  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      libMesh::Number u, v, w, T;
      u = c.interior_value(this->_u_var, qp);
      v = c.interior_value(this->_v_var, qp);
      if (this->_dim == 3)
        w = c.interior_value(this->_w_var, qp);

      T = c.interior_value(this->_T_var, qp);

      libMesh::Gradient grad_u, grad_v, grad_w, grad_T;
      grad_u = c.interior_gradient(this->_u_var, qp);
      grad_v = c.interior_gradient(this->_v_var, qp);
      if (this->_dim == 3)
       grad_w = c.interior_gradient(this->_w_var, qp);

      grad_T = c.interior_gradient( this->_T_var, qp);

      libMesh::NumberVectorValue U(u,v);
      if (this->_dim == 3)
        U(2) = w;

      libMesh::Number divU = grad_u(0) + grad_v(1);
      if (this->_dim == 3)
	divU += grad_w(2);

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) += (-U*grad_T/T + divU)*p_phi[i][qp]*JxW[qp];
	}
    }

  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokes<Mu,SH,TC>::assemble_momentum_time_deriv( bool request_jacobian, 
									 libMesh::FEMContext& c, 
									 libMesh::FEMSystem* system )
{
  // The number of local degrees of freedom in each variable.
  const unsigned int n_u_dofs = c.dof_indices_var[this->_u_var].size();

  // Check number of dofs is same for _u_var, v_var and w_var.
  libmesh_assert (n_u_dofs == c.dof_indices_var[this->_v_var].size());
  if (this->_dim == 3)
    libmesh_assert (n_u_dofs == c.dof_indices_var[this->_w_var].size());

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[this->_u_var]->get_JxW();

  // The pressure shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& u_phi =
    c.element_fe_var[this->_u_var]->get_phi();

// The velocity shape function gradients at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
    c.element_fe_var[this->_u_var]->get_dphi();

  libMesh::DenseSubVector<Number> &Fu = *c.elem_subresiduals[this->_u_var]; // R_{u}
  libMesh::DenseSubVector<Number> &Fv = *c.elem_subresiduals[this->_v_var]; // R_{v}
  libMesh::DenseSubVector<Number> &Fw = *c.elem_subresiduals[this->_w_var]; // R_{w}

  unsigned int n_qpoints = c.element_qrule->n_points();
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      libMesh::Number u, v, w, p, T;
      u = c.interior_value(this->_u_var, qp);
      v = c.interior_value(this->_v_var, qp);
      if (this->_dim == 3)
        w = c.interior_value(this->_w_var, qp);
      p = c.interior_value(this->_p_var, qp);
      T = c.interior_value(this->_T_var, qp);

      libMesh::Gradient grad_u, grad_v, grad_w;
      grad_u = c.interior_gradient(this->_u_var, qp);
      grad_v = c.interior_gradient(this->_v_var, qp);
      if (this->_dim == 3)
       grad_w = c.interior_gradient(this->_w_var, qp);

      libMesh::NumberVectorValue grad_uT( grad_u(0), grad_v(0) ); 
      libMesh::NumberVectorValue grad_vT( grad_u(1), grad_v(1) );
      libMesh::NumberVectorValue grad_wT;
      if( this->_dim == 3 )
	{
	  grad_uT(2) = grad_w(0);
	  grad_vT(2) = grad_w(1);
	  grad_wT = libMesh::NumberVectorValue( grad_u(2), grad_v(2), grad_w(2) );
	}

      libMesh::NumberVectorValue U(u,v);
      if (this->_dim == 3)
        U(2) = w;

      libMesh::Number divU = grad_u(0) + grad_v(1);
      if (this->_dim == 3)
	divU += grad_w(2);

      libMesh::Number rho = this->compute_rho( T, this->get_p0_steady(c,qp) );
      
      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fu(i) += ( -rho*U*grad_u*u_phi[i][qp]                 // convection term
		     + p*u_gradphi[i][qp](0)                           // pressure term
		     - this->_mu(T)*(u_gradphi[i][qp]*grad_u + u_gradphi[i][qp]*grad_uT
			       - 2.0/3.0*divU*u_gradphi[i][qp](0) )    // diffusion term
		     + rho*this->_g(0)*u_phi[i][qp]                 // hydrostatic term
		     )*JxW[qp]; 

          Fv(i) += ( -rho*U*grad_v*u_phi[i][qp]                 // convection term
		     + p*u_gradphi[i][qp](1)                           // pressure term
		     - this->_mu(T)*(u_gradphi[i][qp]*grad_v + u_gradphi[i][qp]*grad_vT
			       - 2.0/3.0*divU*u_gradphi[i][qp](1) )    // diffusion term
		     + rho*this->_g(1)*u_phi[i][qp]                 // hydrostatic term
		     )*JxW[qp];
          if (this->_dim == 3)
            {
              Fw(i) += ( -rho*U*grad_w*u_phi[i][qp]                 // convection term
			 + p*u_gradphi[i][qp](2)                           // pressure term
			 - this->_mu(T)*(u_gradphi[i][qp]*grad_w + u_gradphi[i][qp]*grad_wT
				   - 2.0/3.0*divU*u_gradphi[i][qp](2) )    // diffusion term
			 + rho*this->_g(2)*u_phi[i][qp]                 // hydrostatic term
			 )*JxW[qp];
            }

	  /*
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
	  */
	} // End of DoF loop i
    } // End quadrature loop qp

  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokes<Mu,SH,TC>::assemble_energy_time_deriv( bool request_jacobian, 
								       libMesh::FEMContext& c, 
								       libMesh::FEMSystem* system )
{
  // The number of local degrees of freedom in each variable.
  const unsigned int n_T_dofs = c.dof_indices_var[this->_T_var].size();

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[this->_T_var]->get_JxW();

  // The temperature shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& T_phi =
    c.element_fe_var[this->_T_var]->get_phi();

  // The temperature shape functions gradients at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
    c.element_fe_var[this->_T_var]->get_dphi();

  libMesh::DenseSubVector<Number> &FT = *c.elem_subresiduals[this->_T_var]; // R_{T}

  unsigned int n_qpoints = c.element_qrule->n_points();
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      libMesh::Number u, v, w, T;
      u = c.interior_value(this->_u_var, qp);
      v = c.interior_value(this->_v_var, qp);
      if (this->_dim == 3)
        w = c.interior_value(this->_w_var, qp);
      T = c.interior_value(this->_T_var, qp);

      libMesh::Gradient grad_T;
      grad_T = c.interior_gradient(this->_T_var, qp);

      libMesh::NumberVectorValue U(u,v);
      if (this->_dim == 3)
        U(2) = w;

      libMesh::Number k = this->_k(T);
      libMesh::Number cp = this->_cp(T);

      libMesh::Number rho = this->compute_rho( T, this->get_p0_steady(c,qp) );

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_T_dofs; i++)
        {
          FT(i) += ( -rho*cp*U*grad_T*T_phi[i][qp] // convection term
		     - k*grad_T*T_gradphi[i][qp]            // diffusion term
		     )*JxW[qp]; 
	}
    }

  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokes<Mu,SH,TC>::assemble_continuity_mass_residual( bool request_jacobian, 
									      libMesh::FEMContext& c, 
									      libMesh::FEMSystem* system )
{
  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = 
    c.element_fe_var[this->_u_var]->get_JxW();

  // The shape functions at interior quadrature points.
  const std::vector<std::vector<Real> >& p_phi = 
    c.element_fe_var[this->_p_var]->get_phi();
  
  // The number of local degrees of freedom in each variable
  const unsigned int n_p_dofs = c.dof_indices_var[this->_p_var].size();

  // The subvectors and submatrices we need to fill:
  DenseSubVector<Real> &F_p = *c.elem_subresiduals[this->_p_var];

  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp = 0; qp != n_qpoints; ++qp)
    {
      // For the mass residual, we need to be a little careful.
      // The time integrator is handling the time-discretization
      // for us so we need to supply M(u_fixed)*u for the residual.
      // u_fixed will be given by the fixed_interior_* functions
      // while u will be given by the interior_* functions.
      Real T_dot = c.interior_value(this->_T_var, qp);

      Real T = c.fixed_interior_value(this->_T_var, qp);

      for (unsigned int i = 0; i != n_p_dofs; ++i)
        {
	  F_p(i) += T_dot/T*p_phi[i][qp]*JxW[qp];
	} // End DoF loop i

    } // End quadrature loop qp

  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokes<Mu,SH,TC>::assemble_momentum_mass_residual( bool request_jacobian, 
									    libMesh::FEMContext& c, 
									    libMesh::FEMSystem* system )
{
  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = 
    c.element_fe_var[this->_u_var]->get_JxW();

  // The shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& u_phi = 
    c.element_fe_var[this->_u_var]->get_phi();
  
  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[this->_u_var].size();

  // for convenience
  if (this->_dim != 3)
    this->_w_var = this->_u_var;

  // The subvectors and submatrices we need to fill:
  DenseSubVector<Real> &F_u = *c.elem_subresiduals[this->_u_var];
  DenseSubVector<Real> &F_v = *c.elem_subresiduals[this->_v_var];
  DenseSubVector<Real> &F_w = *c.elem_subresiduals[this->_w_var];

  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp = 0; qp != n_qpoints; ++qp)
    {
      // For the mass residual, we need to be a little careful.
      // The time integrator is handling the time-discretization
      // for us so we need to supply M(u_fixed)*u for the residual.
      // u_fixed will be given by the fixed_interior_* functions
      // while u will be given by the interior_* functions.
      Real u_dot = c.interior_value(this->_u_var, qp);
      Real v_dot = c.interior_value(this->_v_var, qp);

      Real w_dot = 0.0;
      if( this->_dim == 3 )
	Real w_dot = c.interior_value(this->_w_var, qp);

      Real T = c.fixed_interior_value(this->_T_var, qp);
      
      libMesh::Number rho = this->compute_rho(T, this->get_p0_transient(c, qp));
      
      for (unsigned int i = 0; i != n_u_dofs; ++i)
        {
	  F_u(i) += rho*u_dot*u_phi[i][qp]*JxW[qp];
	  F_v(i) += rho*v_dot*u_phi[i][qp]*JxW[qp];

	  if( this->_dim == 3 )
	    F_w(i) += rho*w_dot*u_phi[i][qp]*JxW[qp];
	  
	  /*
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
		  
		} // End DoF loop j
	    } // End Jacobian check
	  */

	} // End DoF loop i
    } // End quadrature loop qp

  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokes<Mu,SH,TC>::assemble_energy_mass_residual( bool request_jacobian, 
									  libMesh::FEMContext& c, 
									  libMesh::FEMSystem* system )
{
  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = 
    c.element_fe_var[this->_T_var]->get_JxW();

  // The shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& T_phi = 
    c.element_fe_var[this->_T_var]->get_phi();
  
  // The number of local degrees of freedom in each variable
  const unsigned int n_T_dofs = c.dof_indices_var[this->_T_var].size();

  // The subvectors and submatrices we need to fill:
  DenseSubVector<Real> &F_T = *c.elem_subresiduals[this->_u_var];

  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp = 0; qp != n_qpoints; ++qp)
    {
      // For the mass residual, we need to be a little careful.
      // The time integrator is handling the time-discretization
      // for us so we need to supply M(u_fixed)*u for the residual.
      // u_fixed will be given by the fixed_interior_* functions
      // while u will be given by the interior_* functions.
      Real T_dot = c.interior_value(this->_T_var, qp);

      Real T = c.fixed_interior_value(this->_T_var, qp);

      Real cp = this->_cp(T);
      
      libMesh::Number rho = this->compute_rho(T, this->get_p0_transient(c, qp));
      
      for (unsigned int i = 0; i != n_T_dofs; ++i)
        {
	  F_T(i) += rho*cp*T_dot*T_phi[i][qp]*JxW[qp];
	} // End DoF loop i

    } // End quadrature loop qp

  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokes<Mu,SH,TC>::assemble_thermo_press_elem_time_deriv( bool request_jacobian, 
										  libMesh::FEMContext& c, 
										  libMesh::FEMSystem* system )
{
  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = 
    c.element_fe_var[this->_T_var]->get_JxW();

  // The number of local degrees of freedom in each variable
  const unsigned int n_p0_dofs = c.dof_indices_var[this->_p0_var].size();

  // The subvectors and submatrices we need to fill:
  DenseSubVector<Real> &F_p0 = *c.elem_subresiduals[this->_p0_var];

  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp = 0; qp != n_qpoints; ++qp)
    {
      libMesh::Number T;
      T = c.interior_value(this->_T_var, qp);

      libMesh::Gradient grad_u, grad_v, grad_w;
      grad_u = c.interior_gradient(this->_u_var, qp);
      grad_v = c.interior_gradient(this->_v_var, qp);
      if (this->_dim == 3)
       grad_w = c.interior_gradient(this->_w_var, qp);

      libMesh::Number divU = grad_u(0) + grad_v(1);
      if(this->_dim==3)
	divU += grad_w(2);

      libMesh::Number cp = this->_cp(T);
      libMesh::Number cv = cp + this->_R;
      libMesh::Number gamma = cp/cv;
      libMesh::Number gamma_ratio = gamma/(gamma-1.0);

      libMesh::Number p0 = c.interior_value( this->_p0_var, qp );

      for (unsigned int i = 0; i != n_p0_dofs; ++i)
        {
	  F_p0(i) += (p0/T - this->_p0/this->_T0)*JxW[qp];
	  //F_p0(i) -= p0*gamma_ratio*divU*JxW[qp];
	} // End DoF loop i
    }

  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokes<Mu,SH,TC>::assemble_thermo_press_side_time_deriv( bool request_jacobian, 
										  libMesh::FEMContext& c, 
										  libMesh::FEMSystem* system )
{
  // The number of local degrees of freedom in each variable.
  const unsigned int n_p0_dofs = c.dof_indices_var[this->_p0_var].size();

  // Element Jacobian * quadrature weight for side integration.
  const std::vector<libMesh::Real> &JxW_side = c.side_fe_var[this->_T_var]->get_JxW();

  const std::vector<Point> &normals = c.side_fe_var[this->_T_var]->get_normals();

  libMesh::DenseSubVector<Number> &F_p0 = *c.elem_subresiduals[this->_p0_var]; // residual

  // Physical location of the quadrature points
  const std::vector<libMesh::Point>& qpoint =
    c.side_fe_var[this->_T_var]->get_xyz();

  unsigned int n_qpoints = c.side_qrule->n_points();
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      libMesh::Number T = c.side_value( this->_T_var, qp );
      libMesh::Gradient U = ( c.side_value( this->_u_var, qp ),
			      c.side_value( this->_v_var, qp ) );
      libMesh::Gradient grad_T = c.side_gradient( this->_T_var, qp );

      libMesh::Number p0 = c.side_value( this->_p0_var, qp );

      libMesh::Number k = this->_k(T);
      libMesh::Number cp = this->_cp(T);

      libMesh::Number cv = cp + this->_R;
      libMesh::Number gamma = cp/cv;
      libMesh::Number gamma_ratio = gamma/(gamma-1.0);

      //std::cout << "U = " << U << std::endl;

      //std::cout << "x = " << qpoint[qp] << ", grad_T = " << grad_T << std::endl;

      for (unsigned int i=0; i != n_p0_dofs; i++)
	{
	  //F_p0(i) += (k*grad_T*normals[qp] - p0*gamma_ratio*U*normals[qp]  )*JxW_side[qp];
	}
    }
  
  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokes<Mu,SH,TC>::assemble_thermo_press_mass_residual( bool request_jacobian, 
										libMesh::FEMContext& c, 
										libMesh::FEMSystem* system )
{
  // The number of local degrees of freedom in each variable.
  const unsigned int n_p0_dofs = c.dof_indices_var[this->_p0_var].size();
  const unsigned int n_T_dofs = c.dof_indices_var[this->_T_var].size();
  const unsigned int n_p_dofs = c.dof_indices_var[this->_p_var].size();

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = 
    c.element_fe_var[this->_T_var]->get_JxW();

  // The temperature shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& T_phi =
    c.element_fe_var[this->_T_var]->get_phi();

  // The temperature shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& p_phi =
    c.element_fe_var[this->_p_var]->get_phi();

  // The subvectors and submatrices we need to fill:
  DenseSubVector<Real> &F_p0 = *c.elem_subresiduals[this->_p0_var];
  DenseSubVector<Real> &F_T = *c.elem_subresiduals[this->_T_var];
  DenseSubVector<Real> &F_p = *c.elem_subresiduals[this->_p_var];

  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp = 0; qp != n_qpoints; ++qp)
    {
      libMesh::Number T;
      T = c.fixed_interior_value(this->_T_var, qp);

      libMesh::Number cp = this->_cp(T);
      libMesh::Number cv = cp + this->_R;
      libMesh::Number gamma = cp/cv;
      libMesh::Number one_over_gamma = 1.0/(gamma-1.0);

      libMesh::Number p0_dot = c.interior_value(this->_p0_var, qp );

      libMesh::Number p0 = c.fixed_interior_value(this->_p0_var, qp );

      for (unsigned int i=0; i != n_p0_dofs; i++)
	{
	  F_p0(i) += p0_dot*one_over_gamma*JxW[qp];
	}

      for (unsigned int i=0; i != n_T_dofs; i++)
	{
	  F_T(i) -= p0_dot*T_phi[i][qp]*JxW[qp];
	}

      for (unsigned int i=0; i != n_p_dofs; i++)
	{
	  F_p(i) -= p0_dot/p0*p_phi[i][qp]*JxW[qp];
	}

    }
  return;
}

// Instantiate
template class GRINS::LowMachNavierStokes<GRINS::ConstantViscosity,GRINS::ConstantSpecificHeat,GRINS::ConstantConductivity>;
