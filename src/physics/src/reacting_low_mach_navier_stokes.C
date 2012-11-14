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

#include "reacting_low_mach_navier_stokes.h"

namespace GRINS
{
  template<class Mixture>
  ReactingLowMachNavierStokes<Mixture>::ReactingLowMachNavierStokes(const PhysicsName& physics_name, const GetPot& input)
    : ReactingLowMachNavierStokesBase<Mixture>(physics_name,input),
      _p_pinning(input,physics_name)
  {
    this->read_input_options(input);

    return;
  }

  template<class Mixture>
  ReactingLowMachNavierStokes<Mixture>::~ReactingLowMachNavierStokes()
  {
    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::read_input_options( const GetPot& input )
  {
    // Other quantities read in base class

    // Read pressure pinning information
    this->_pin_pressure = input("Physics/"+reacting_low_mach_navier_stokes+"/pin_pressure", false );
  
    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::init_context( libMesh::DiffContext &context )
  {
    // First call base class
    GRINS::ReactingLowMachNavierStokesBase<Mixture>::init_context(context);

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

  template<class Mixture>
  bool ReactingLowMachNavierStokes<Mixture>::element_time_derivative( bool request_jacobian,
								      libMesh::DiffContext& context,
								      libMesh::FEMSystem* system )
  {
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);
    
    unsigned int n_qpoints = c.element_qrule->n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	// Compute species mass fractions at quadrature points
	std::vector<Real> Y;
	Y.reserve(this->_n_species);
	for( unsigned int s = 0; s < this->_n_species; s++ )
	  {
	    Y.push_back( c.interior_value(this->_species_vars[s],qp) );
	  }

	// Build up cache at element interior quadrature point
	ReactingFlowCache cache( c.interior_value(this->_T_var, qp), 
				 this->get_p0_steady(c,qp), Y );

	this->build_reacting_flow_cache(c, cache, qp);

	this->assemble_mass_time_deriv(c, cache, qp);
	this->assemble_species_time_deriv(c, cache, qp);
	this->assemble_momentum_time_deriv(c, cache, qp);
	this->assemble_energy_time_deriv(c, cache, qp);
      }

    // Pin p = p_value at p_point
    if( this->_pin_pressure )
      {
	this->_p_pinning.pin_value( context, request_jacobian, this->_p_var );
      }

    return request_jacobian;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::assemble_mass_time_deriv(libMesh::FEMContext& c, 
								      const ReactingFlowCache& cache, 
								      unsigned int qp)
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
    
    libMesh::Number T = cache.T();
    const libMesh::NumberVectorValue& U = cache.U();
    libMesh::Number divU = cache.divU();
    const libMesh::Gradient& grad_T = cache.grad_T();
    
    for (unsigned int i=0; i != n_p_dofs; i++)
      {
	Fp(i) += MOLAR_MASS_TERM (-U*grad_T/T + divU)*p_phi[i][qp]*JxW[qp];
      }

    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::assemble_species_time_deriv(libMesh::FEMContext& c, 
									 const ReactingFlowCache& cache, 
									 unsigned int qp)
  {
    
    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::assemble_momentum_time_deriv(libMesh::FEMContext& c, 
									  const ReactingFlowCache& cache, 
									  unsigned int qp)
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
    
    libMesh::Number rho = cache.rho();
    const libMesh::NumberVectorValue& U = cache.U();
    libMesh::Number divU = cache.divU();
    libMesh::Number mu = cache.mu();
    libMesh::Number p = cache.p_hydro();
    const libMesh::Gradient& grad_u = cache.grad_u();
    const libMesh::Gradient& grad_v = cache.grad_v();
    const libMesh::Gradient& grad_w = cache.grad_w();

    libMesh::NumberVectorValue grad_uT( grad_u(0), grad_v(0) ); 
    libMesh::NumberVectorValue grad_vT( grad_u(1), grad_v(1) );
    libMesh::NumberVectorValue grad_wT;
    if( this->_dim == 3 )
      {
	grad_uT(2) = grad_w(0);
	grad_vT(2) = grad_w(1);
	grad_wT = libMesh::NumberVectorValue( grad_u(2), grad_v(2), grad_w(2) );
      }

    // Now a loop over the pressure degrees of freedom.  This
    // computes the contributions of the continuity equation.
    for (unsigned int i=0; i != n_u_dofs; i++)
      {
	Fu(i) += ( -rho*U*grad_u*u_phi[i][qp]                 // convection term
		   + p*u_gradphi[i][qp](0)                           // pressure term
		   - mu*(u_gradphi[i][qp]*grad_u + u_gradphi[i][qp]*grad_uT
			 - 2.0/3.0*divU*u_gradphi[i][qp](0) )    // diffusion term
		   + rho*this->_g(0)*u_phi[i][qp]                 // hydrostatic term
		   )*JxW[qp]; 

	Fv(i) += ( -rho*U*grad_v*u_phi[i][qp]                 // convection term
		   + p*u_gradphi[i][qp](1)                           // pressure term
		   - mu*(u_gradphi[i][qp]*grad_v + u_gradphi[i][qp]*grad_vT
			 - 2.0/3.0*divU*u_gradphi[i][qp](1) )    // diffusion term
		   + rho*this->_g(1)*u_phi[i][qp]                 // hydrostatic term
		   )*JxW[qp];

	if (this->_dim == 3)
	  {
	    Fw(i) += ( -rho*U*grad_w*u_phi[i][qp]                 // convection term
		       + p*u_gradphi[i][qp](2)                           // pressure term
		       - mu*(u_gradphi[i][qp]*grad_w + u_gradphi[i][qp]*grad_wT
			     - 2.0/3.0*divU*u_gradphi[i][qp](2) )    // diffusion term
		       + rho*this->_g(2)*u_phi[i][qp]                 // hydrostatic term
		       )*JxW[qp];
	  }
      }
    return;
  }

  template<class Mixture>
  void ReactingLowMachNavierStokes<Mixture>::assemble_energy_time_deriv(libMesh::FEMContext& c, 
									const ReactingFlowCache& cache, 
									unsigned int qp)
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

    libMesh::Number rho = cache.rho() ;
    const libMesh::NumberVectorValue& U = cache.U();
    libMesh::Number cp = cache.cp();
    libMesh::Number k = cache.k();
    const libMesh::Gradient& grad_T = cache.grad_T();
    
    for (unsigned int i=0; i != n_T_dofs; i++)
      {
	FT(i) += ( -rho*cp*U*grad_T*T_phi[i][qp] // convection term
		   - k*grad_T*T_gradphi[i][qp]            // diffusion term
		   )*JxW[qp]; 
      }

    return;
  }

  // Instantiate
#ifdef HAVE_CANTERA
  template class ReactingLowMachNavierStokes< IdealGasMixture<CanteraThermodynamics,CanteraTransport,CanteraKinetics> >;
#endif

} // namespace GRINS
