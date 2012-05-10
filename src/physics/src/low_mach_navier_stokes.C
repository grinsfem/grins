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
  : Physics(physics_name),
    _p_pinning(input,physics_name)
{
  this->read_input_options(input);

  // This is deleted in the base class
  _bc_handler = new GRINS::LowMachNavierStokesBCHandling( physics_name, input );

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
  // Read FE info
  this->_V_FE_family =
    libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+low_mach_navier_stokes+"/V_FE_family", "LAGRANGE") );

  this->_P_FE_family =
    libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+low_mach_navier_stokes+"/P_FE_family", "LAGRANGE") );

  this->_T_FE_family =
    libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+low_mach_navier_stokes+"/T_FE_family", "LAGRANGE") );

  this->_V_order =
    libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+low_mach_navier_stokes+"/V_order", "SECOND") );

  this->_P_order =
    libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+low_mach_navier_stokes+"/P_order", "FIRST") );

  this->_T_order =
    libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+low_mach_navier_stokes+"/T_order", "SECOND") );

  // Read variable naming info
  this->_u_var_name = input("Physics/VariableNames/u_velocity", GRINS::u_var_name_default );
  this->_v_var_name = input("Physics/VariableNames/v_velocity", GRINS::v_var_name_default );
  this->_w_var_name = input("Physics/VariableNames/w_velocity", GRINS::w_var_name_default );
  this->_p_var_name = input("Physics/VariableNames/pressure", GRINS::p_var_name_default );
  this->_T_var_name = input("Physics/VariableNames/temperature", GRINS::T_var_name_default );

  // Read material parameters
  this->_mu.read_input_options( input );
  this->_cp.read_input_options( input );
  this->_k.read_input_options( input );

  // Read thermodynamic state info
  _p0 = input("Physics/"+low_mach_navier_stokes+"/p0", 0.0 ); /* thermodynamic pressure */
  _R  = input("Physics/"+low_mach_navier_stokes+"/R", 0.0 ); /* gas constant */

  if( _R <= 0.0 )
    {
      std::cerr << "=========================================" << std::endl
		<< " Error: Gas constant R must be positive. " << std::endl
		<< " Detected value R = " << _R << std::endl
		<< "=========================================" << std::endl;
      libmesh_error();
    }

  _p0_over_R = _p0/_R;

  _enable_thermo_press_calc = input("Physics/"+low_mach_navier_stokes+"/enable_thermo_press_calc", false );

  if( _enable_thermo_press_calc )
    {
      _p0_var_name = input("Physics/VariableNames/thermo_presure", "P0" );
    }

  // Read gravity vector
  unsigned int g_dim = input.vector_variable_size("Physics/"+low_mach_navier_stokes+"/g");

  _g(0) = input("Physics/"+low_mach_navier_stokes+"/g", 0.0, 0 );
  _g(1) = input("Physics/"+low_mach_navier_stokes+"/g", 0.0, 1 );
  
  if( g_dim == 3)
    _g(2) = input("Physics/"+low_mach_navier_stokes+"/g", 0.0, 2 );

  // Read pressure pinning information
  _pin_pressure = input("Physics/"+low_mach_navier_stokes+"/pin_pressure", true );
  
  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokes<Mu,SH,TC>::init_variables( libMesh::FEMSystem* system )
{
  // Get libMesh to assign an index for each variable
  this->_dim = system->get_mesh().mesh_dimension();

  _u_var = system->add_variable( _u_var_name, this->_V_order, _V_FE_family);
  _v_var = system->add_variable( _v_var_name, this->_V_order, _V_FE_family);

  if (_dim == 3)
    _w_var = system->add_variable( _w_var_name, this->_V_order, _V_FE_family);

  _p_var = system->add_variable( _p_var_name, this->_P_order, _P_FE_family);
  _T_var = system->add_variable( _T_var_name, this->_T_order, _T_FE_family);

  /* If we need to compute the thermodynamic pressure, we force this to be a first
     order scalar variable. */
  if( _enable_thermo_press_calc )
    _p0_var = system->add_variable( _p0_var_name, FIRST, SCALAR);

  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokes<Mu,SH,TC>::set_time_evolving_vars( libMesh::FEMSystem* system )
{
  const unsigned int dim = system->get_mesh().mesh_dimension();

  system->time_evolving(_u_var);
  system->time_evolving(_v_var);

  if (dim == 3)
    system->time_evolving(_w_var);

  system->time_evolving(_p_var);
  system->time_evolving(_T_var);

  return;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokes<Mu,SH,TC>::init_context( libMesh::DiffContext &context )
{
  libMesh::FEMContext &c = libmesh_cast_ref<libMesh::FEMContext&>(context);

  // We should prerequest all the data
  // we will need to build the linear system
  // or evaluate a quantity of interest.
  c.element_fe_var[_u_var]->get_JxW();
  c.element_fe_var[_u_var]->get_phi();
  c.element_fe_var[_u_var]->get_dphi();
  c.element_fe_var[_u_var]->get_xyz();

  c.element_fe_var[_T_var]->get_JxW();
  c.element_fe_var[_T_var]->get_phi();
  c.element_fe_var[_T_var]->get_dphi();
  c.element_fe_var[_T_var]->get_xyz();

  c.element_fe_var[_p_var]->get_phi();
  c.element_fe_var[_p_var]->get_xyz();

  c.side_fe_var[_u_var]->get_JxW();
  c.side_fe_var[_u_var]->get_phi();
  c.side_fe_var[_u_var]->get_dphi();
  c.side_fe_var[_u_var]->get_xyz();

  c.side_fe_var[_T_var]->get_JxW();
  c.side_fe_var[_T_var]->get_phi();
  c.side_fe_var[_T_var]->get_dphi();
  c.side_fe_var[_T_var]->get_xyz();

  return;
}

template<class Mu, class SH, class TC>
bool GRINS::LowMachNavierStokes<Mu,SH,TC>::element_time_derivative( bool request_jacobian,
								    libMesh::DiffContext& context,
								    libMesh::FEMSystem* system )
{
  // Testing phase. We'll rely on numerical Jacobians.
  request_jacobian = false;

#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("LowMachNavierStokes::element_time_derivative");
#endif

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  this->assemble_mass_time_deriv( request_jacobian, c, system );
  this->assemble_momentum_time_deriv( request_jacobian, c, system );
  this->assemble_energy_time_deriv( request_jacobian, c, system );

  // Pin p = p_value at p_point
  if( _pin_pressure )
    {
      _p_pinning.pin_value( context, request_jacobian, _p_var);
    }

  if( _enable_thermo_press_calc )
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
  if( _enable_thermo_press_calc )
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
  // Testing phase. We'll rely on numerical Jacobians.
  request_jacobian = false;

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  this->assemble_continuity_mass_residual( request_jacobian, c, system );

  this->assemble_momentum_mass_residual( request_jacobian, c, system );

  this->assemble_energy_mass_residual( request_jacobian, c, system );

  if( _enable_thermo_press_calc )
    this->assemble_thermo_press_mass_residual( request_jacobian, c, system );

  return request_jacobian;
}

template<class Mu, class SH, class TC>
void GRINS::LowMachNavierStokes<Mu,SH,TC>::assemble_mass_time_deriv( bool request_jacobian, 
								     libMesh::FEMContext& c, 
								     libMesh::FEMSystem* system )
{
  // The number of local degrees of freedom in each variable.
  const unsigned int n_p_dofs = c.dof_indices_var[_p_var].size();

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[_u_var]->get_JxW();

  // The pressure shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& p_phi =
    c.element_fe_var[_p_var]->get_phi();

  libMesh::DenseSubVector<Number> &Fp = *c.elem_subresiduals[_p_var]; // R_{p}

  unsigned int n_qpoints = c.element_qrule->n_points();
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      libMesh::Number u, v, w, T;
      u = c.interior_value(_u_var, qp);
      v = c.interior_value(_v_var, qp);
      if (_dim == 3)
        w = c.interior_value(_w_var, qp);

      T = c.interior_value(_T_var, qp);

      libMesh::Gradient grad_u, grad_v, grad_w, grad_T;
      grad_u = c.interior_gradient(_u_var, qp);
      grad_v = c.interior_gradient(_v_var, qp);
      if (_dim == 3)
       grad_w = c.interior_gradient(_w_var, qp);

      grad_T = c.interior_gradient( _T_var, qp);

      libMesh::NumberVectorValue U(u,v);
      if (_dim == 3)
        U(2) = w;

      libMesh::Number divU = grad_u(0) + grad_v(1);
      if (_dim == 3)
	divU += grad_w(2);

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) += (-U*grad_T + T*divU)*p_phi[i][qp]*JxW[qp];
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
  const unsigned int n_u_dofs = c.dof_indices_var[_u_var].size();

  // Check number of dofs is same for _u_var, v_var and w_var.
  libmesh_assert (n_u_dofs == c.dof_indices_var[_v_var].size());
  if (_dim == 3)
    libmesh_assert (n_u_dofs == c.dof_indices_var[_w_var].size());

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[_u_var]->get_JxW();

  // The pressure shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& u_phi =
    c.element_fe_var[_u_var]->get_phi();

// The velocity shape function gradients at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
    c.element_fe_var[_u_var]->get_dphi();

  libMesh::DenseSubVector<Number> &Fu = *c.elem_subresiduals[_u_var]; // R_{u}
  libMesh::DenseSubVector<Number> &Fv = *c.elem_subresiduals[_v_var]; // R_{v}
  libMesh::DenseSubVector<Number> &Fw = *c.elem_subresiduals[_w_var]; // R_{w}

  unsigned int n_qpoints = c.element_qrule->n_points();
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      libMesh::Number u, v, w, p, T;
      u = c.interior_value(_u_var, qp);
      v = c.interior_value(_v_var, qp);
      if (_dim == 3)
        w = c.interior_value(_w_var, qp);
      p = c.interior_value(_p_var, qp);
      T = c.interior_value(_T_var, qp);

      libMesh::Gradient grad_u, grad_v, grad_w;
      grad_u = c.interior_gradient(_u_var, qp);
      grad_v = c.interior_gradient(_v_var, qp);
      if (_dim == 3)
       grad_w = c.interior_gradient(_w_var, qp);

      libMesh::NumberVectorValue U(u,v);
      if (_dim == 3)
        U(2) = w;

      libMesh::Number divU = grad_u(0) + grad_v(1);
      if (_dim == 3)
	divU += grad_w(2);

      libMesh::Number p0_over_R;
      if( _enable_thermo_press_calc )
	{
	  libMesh::Number p0 = c.fixed_interior_value( _p0_var,qp );
	  p0_over_R = p0/_R;
	}
      else
	{
	  p0_over_R = _p0_over_R;
	}

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fu(i) += ( -p0_over_R*U*grad_u*u_phi[i][qp]                 // convection term
		     + p*u_gradphi[i][qp](0)                           // pressure term
		     - _mu(T)*(grad_u*u_gradphi[i][qp] 
			       - 2.0/3.0*divU*u_gradphi[i][qp](0) )    // diffusion term
		     + p0_over_R/T*_g(0)*u_phi[i][qp]                 // hydrostatic term
		     )*JxW[qp]; 

          Fv(i) += ( -p0_over_R*U*grad_v*u_phi[i][qp]                 // convection term
		     + p*u_gradphi[i][qp](1)                           // pressure term
		     - _mu(T)*(grad_v*u_gradphi[i][qp] 
			       - 2.0/3.0*divU*u_gradphi[i][qp](1) )    // diffusion term
		     + p0_over_R/T*_g(1)*u_phi[i][qp]                 // hydrostatic term
		     )*JxW[qp];
          if (_dim == 3)
            {
              Fw(i) += ( -p0_over_R*U*grad_w*u_phi[i][qp]                 // convection term
			 + p*u_gradphi[i][qp](2)                           // pressure term
			 - _mu(T)*(grad_w*u_gradphi[i][qp] 
				   - 2.0/3.0*divU*u_gradphi[i][qp](2) )    // diffusion term
			 + p0_over_R/T*_g(2)*u_phi[i][qp]                 // hydrostatic term
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
  const unsigned int n_T_dofs = c.dof_indices_var[_T_var].size();

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[_T_var]->get_JxW();

  // The temperature shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& T_phi =
    c.element_fe_var[_T_var]->get_phi();

  // The temperature shape functions gradients at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
    c.element_fe_var[_T_var]->get_dphi();

  libMesh::DenseSubVector<Number> &FT = *c.elem_subresiduals[_T_var]; // R_{T}

  unsigned int n_qpoints = c.element_qrule->n_points();
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      libMesh::Number u, v, w, T;
      u = c.interior_value(_u_var, qp);
      v = c.interior_value(_v_var, qp);
      if (_dim == 3)
        w = c.interior_value(_w_var, qp);
      T = c.interior_value(_T_var, qp);

      libMesh::Gradient grad_T;
      grad_T = c.interior_gradient(_T_var, qp);

      libMesh::NumberVectorValue U(u,v);
      if (_dim == 3)
        U(2) = w;

      libMesh::Number k = this->_k(T);
      libMesh::Number cp = this->_cp(T);

      libMesh::Number p0_over_R;
      if( _enable_thermo_press_calc )
	{
	  libMesh::Number p0 = c.fixed_interior_value( _p0_var,qp );
	  p0_over_R = p0/_R;
	}
      else
	{
	  p0_over_R = _p0_over_R;
	}

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_T_dofs; i++)
        {
          FT(i) += ( -p0_over_R*cp/T*U*grad_T*T_phi[i][qp] // convection term
		     + k*grad_T*T_gradphi[i][qp]            // diffusion term
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
    c.element_fe_var[_u_var]->get_JxW();

  // The shape functions at interior quadrature points.
  const std::vector<std::vector<Real> >& p_phi = 
    c.element_fe_var[_p_var]->get_phi();
  
  // The number of local degrees of freedom in each variable
  const unsigned int n_p_dofs = c.dof_indices_var[_p_var].size();

  // The subvectors and submatrices we need to fill:
  DenseSubVector<Real> &F_p = *c.elem_subresiduals[_p_var];

  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp = 0; qp != n_qpoints; ++qp)
    {
      // For the mass residual, we need to be a little careful.
      // The time integrator is handling the time-discretization
      // for us so we need to supply M(u_fixed)*u for the residual.
      // u_fixed will be given by the fixed_interior_* functions
      // while u will be given by the interior_* functions.
      Real T_dot = c.interior_value(_T_var, qp);

      for (unsigned int i = 0; i != n_p_dofs; ++i)
        {
	  F_p(i) += T_dot*p_phi[i][qp]*JxW[qp];
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
    c.element_fe_var[_u_var]->get_JxW();

  // The shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& u_phi = 
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

      Real T = c.fixed_interior_value(_T_var, qp);
      
      libMesh::Number p0_over_R;
      if( _enable_thermo_press_calc )
	{
	  libMesh::Number p0 = c.fixed_interior_value( _p0_var,qp );
	  p0_over_R = p0/_R;
	}
      else
	{
	  p0_over_R = _p0_over_R;
	}

      for (unsigned int i = 0; i != n_u_dofs; ++i)
        {
	  F_u(i) += p0_over_R/T*u_dot*u_phi[i][qp]*JxW[qp];
	  F_v(i) += p0_over_R/T*v_dot*u_phi[i][qp]*JxW[qp];

	  if( _dim == 3 )
	    F_w(i) += p0_over_R/T*w_dot*u_phi[i][qp]*JxW[qp];
	  
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
    c.element_fe_var[_T_var]->get_JxW();

  // The shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& T_phi = 
    c.element_fe_var[_T_var]->get_phi();
  
  // The number of local degrees of freedom in each variable
  const unsigned int n_T_dofs = c.dof_indices_var[_T_var].size();

  // The subvectors and submatrices we need to fill:
  DenseSubVector<Real> &F_T = *c.elem_subresiduals[_u_var];

  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp = 0; qp != n_qpoints; ++qp)
    {
      // For the mass residual, we need to be a little careful.
      // The time integrator is handling the time-discretization
      // for us so we need to supply M(u_fixed)*u for the residual.
      // u_fixed will be given by the fixed_interior_* functions
      // while u will be given by the interior_* functions.
      Real T_dot = c.interior_value(_T_var, qp);

      Real T = c.fixed_interior_value(_T_var, qp);

      Real cp = this->_cp(T);
      
      libMesh::Number p0_over_R;
      if( _enable_thermo_press_calc )
	{
	  libMesh::Number p0 = c.fixed_interior_value( _p0_var,qp );
	  p0_over_R = p0/_R;
	}
      else
	{
	  p0_over_R = _p0_over_R;
	}

      for (unsigned int i = 0; i != n_T_dofs; ++i)
        {
	  F_T(i) += p0_over_R/T*cp*T_dot*T_phi[i][qp]*JxW[qp];
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
    c.element_fe_var[_T_var]->get_JxW();
  
  // The number of local degrees of freedom in each variable
  const unsigned int n_p0_dofs = c.dof_indices_var[_p0_var].size();

  // The subvectors and submatrices we need to fill:
  DenseSubVector<Real> &F_p = *c.elem_subresiduals[_p0_var];

  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp = 0; qp != n_qpoints; ++qp)
    {
      libMesh::Number T;
      T = c.interior_value(_T_var, qp);

      libMesh::Gradient grad_u, grad_v, grad_w;
      grad_u = c.interior_gradient(_u_var, qp);
      grad_v = c.interior_gradient(_v_var, qp);
      if (_dim == 3)
       grad_w = c.interior_gradient(_w_var, qp);

      libMesh::Number divU = grad_u(0) + grad_v(1);
      if(_dim==3)
	divU += grad_w(2);

      libMesh::Number cp = _cp(T);
      libMesh::Number cv = cp + _R;
      libMesh::Number gamma = cp/cv;
      libMesh::Number gamma_ratio = gamma/(gamma-1.0);

      libMesh::Number p0 = c.interior_value( _p0_var, qp );

      for (unsigned int i = 0; i != n_p0_dofs; ++i)
        {
	  F_p(i) -= p0*gamma_ratio*divU*JxW[qp];
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
  const unsigned int n_p0_dofs = c.dof_indices_var[_p0_var].size();

  // Element Jacobian * quadrature weight for side integration.
  const std::vector<libMesh::Real> &JxW_side = c.side_fe_var[_T_var]->get_JxW();

  const std::vector<Point> &normals = c.side_fe_var[_T_var]->get_normals();

  libMesh::DenseSubVector<Number> &F_p = *c.elem_subresiduals[_p0_var]; // residual

  unsigned int n_qpoints = c.side_qrule->n_points();
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      libMesh::Number T = c.side_value( _T_var, qp );
      libMesh::Gradient grad_T = c.side_gradient( _T_var, qp );

      libMesh::Number k = _k(T);

      for (unsigned int i=0; i != n_p0_dofs; i++)
	{
	  F_p(i) -= k*grad_T*normals[qp]*JxW_side[qp];
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
  const unsigned int n_p0_dofs = c.dof_indices_var[_p0_var].size();
  const unsigned int n_T_dofs = c.dof_indices_var[_T_var].size();

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = 
    c.element_fe_var[_T_var]->get_JxW();

  // The temperature shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& T_phi =
    c.element_fe_var[_T_var]->get_phi();

  // The subvectors and submatrices we need to fill:
  DenseSubVector<Real> &F_p = *c.elem_subresiduals[_p0_var];
  DenseSubVector<Real> &F_T = *c.elem_subresiduals[_T_var];

  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp = 0; qp != n_qpoints; ++qp)
    {
      libMesh::Number T;
      T = c.fixed_interior_value(_T_var, qp);

      libMesh::Number cp = _cp(T);
      libMesh::Number cv = cp + _R;
      libMesh::Number gamma = cp/cv;
      libMesh::Number one_over_gamma = 1.0/(gamma-1.0);

      libMesh::Number p0_dot = c.interior_value(_p0_var, qp );

      for (unsigned int i=0; i != n_p0_dofs; i++)
	{
	  F_p(i) += p0_dot*one_over_gamma*JxW[qp];
	}

      for (unsigned int i=0; i != n_T_dofs; i++)
	{
	  F_T(i) -= p0_dot*one_over_gamma*T_phi[i][qp]*JxW[qp];
	}

    }
  return;
}

// Instantiate
template class GRINS::LowMachNavierStokes<GRINS::ConstantViscosity,GRINS::ConstantSpecificHeat,GRINS::ConstantConductivity>;
