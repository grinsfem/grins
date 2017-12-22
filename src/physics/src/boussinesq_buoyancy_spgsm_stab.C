//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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

// This class
#include "grins_config.h"
#include "grins/boussinesq_buoyancy_spgsm_stab.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/inc_nav_stokes_macro.h"
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace GRINS
{
  template<class Mu>
  BoussinesqBuoyancySPGSMStabilization<Mu>::BoussinesqBuoyancySPGSMStabilization( const std::string& physics_name, const GetPot& input )
    : BoussinesqBuoyancyBase(physics_name,input),
      _flow_stab_helper(physics_name+"FlowStabHelper", input),
      _temp_stab_helper(physics_name+"TempStabHelper", input),
      _Cp(0.0),
      _k(1.0),
      _mu(input,MaterialsParsing::material_name(input,PhysicsNaming::boussinesq_buoyancy()))
  {
    MaterialsParsing::read_property( input,
                                     "SpecificHeat",
                                     PhysicsNaming::boussinesq_buoyancy(),
                                     (*this),
                                     this->_Cp );

    this->set_parameter
      (_k, input, "Physics/"+PhysicsNaming::heat_transfer()+"/k", _k);
  }

  template<class Mu>
  BoussinesqBuoyancySPGSMStabilization<Mu>::~BoussinesqBuoyancySPGSMStabilization()
  {
    return;
  }

  template<class Mu>
  void BoussinesqBuoyancySPGSMStabilization<Mu>::element_time_derivative
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.get_dof_indices(_flow_vars.u()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(_flow_vars.u())->get_JxW();

    /*
      const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();
    */

    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_flow_vars.u())->get_dphi();

    // Get residuals
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(_flow_vars.u()); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(_flow_vars.v()); // R_{v}
    libMesh::DenseSubVector<libMesh::Number> *Fw = NULL;

    if(this->_flow_vars.dim() == 3)
      {
        Fw = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}
      }

    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.get_element_qrule().n_points();

    libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u());

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::RealGradient g = this->_flow_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_flow_stab_helper.compute_G( fe, context, qp );

        libMesh::RealGradient U( context.interior_value( this->_flow_vars.u(), qp ),
                                 context.interior_value( this->_flow_vars.v(), qp ) );
        if( this->_flow_vars.dim() == 3 )
          {
            U(2) = context.interior_value( this->_flow_vars.w(), qp );
          }

        // Compute the viscosity at this qp
        libMesh::Real mu_qp = this->_mu(context, qp);

        libMesh::Real tau_M = this->_flow_stab_helper.compute_tau_momentum( context, qp, g, G, this->_rho, U, mu_qp, this->_is_steady );

        //libMesh::Real tau_E = this->_temp_stab_helper.compute_tau_energy( context, G, _rho, _Cp, _k,  U, this->_is_steady );

        //libMesh::Real RE = this->_temp_stab_helper.compute_res_energy_steady( context, qp, _rho, _Cp, _k );

        // Compute the solution & its gradient at the old Newton iterate.
        libMesh::Number T;
        T = context.interior_value(_temp_vars.T(), qp);

        libMesh::RealGradient residual = _rho*_beta_T*(T-_T_ref)*_g;

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) += ( -tau_M*residual(0)*_rho*U*u_gradphi[i][qp] )*JxW[qp];
            // + _rho*_beta_T*tau_E*RE*_g(0)*u_phi[i][qp] )*JxW[qp];

            Fv(i) += ( -tau_M*residual(1)*_rho*U*u_gradphi[i][qp] )*JxW[qp];
            // + _rho*_beta_T*tau_E*RE*_g(1)*u_phi[i][qp] )*JxW[qp];

            if (this->_flow_vars.dim() == 3)
              {
                (*Fw)(i) += ( -tau_M*residual(2)*_rho*U*u_gradphi[i][qp] )*JxW[qp];
                // + _rho*_beta_T*tau_E*RE*_g(2)*u_phi[i][qp] )*JxW[qp];
              }

            if (compute_jacobian)
              {
                libmesh_not_implemented();
              } // End compute_jacobian check

          } // End i dof loop
      } // End quadrature loop
  }

  template<class Mu>
  void BoussinesqBuoyancySPGSMStabilization<Mu>::element_constraint
  ( bool compute_jacobian,
    AssemblyContext & context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.get_dof_indices(_press_var.p()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(_flow_vars.u())->get_JxW();

    const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
      context.get_element_fe(this->_press_var.p())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(this->_press_var.p()); // R_{p}

    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.get_element_qrule().n_points();

    libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u());

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::RealGradient g = this->_flow_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_flow_stab_helper.compute_G( fe, context, qp );

        libMesh::RealGradient U( context.interior_value( this->_flow_vars.u(), qp ),
                                 context.interior_value( this->_flow_vars.v(), qp ) );
        if( this->_flow_vars.dim() == 3 )
          U(2) = context.interior_value( this->_flow_vars.w(), qp );

        // Compute the viscosity at this qp
        libMesh::Real mu_qp = this->_mu(context, qp);

        libMesh::Real tau_M = this->_flow_stab_helper.compute_tau_momentum( context, qp, g, G, this->_rho, U, mu_qp, this->_is_steady );

        // Compute the solution & its gradient at the old Newton iterate.
        libMesh::Number T;
        T = context.interior_value(_temp_vars.T(), qp);

        libMesh::RealGradient residual = _rho*_beta_T*(T-_T_ref)*_g;

        // First, an i-loop over the velocity degrees of freedom.
        // We know that n_u_dofs == n_v_dofs so we can compute contributions
        // for both at the same time.
        for (unsigned int i=0; i != n_p_dofs; i++)
          {
            Fp(i) += -tau_M*residual*p_dphi[i][qp]*JxW[qp];

            if (compute_jacobian)
              {
                libmesh_not_implemented();
              } // End compute_jacobian check

          } // End i dof loop

      } // End quadrature loop
  }

  template<class Mu>
  void BoussinesqBuoyancySPGSMStabilization<Mu>::mass_residual( bool /*compute_jacobian*/,
                                                                AssemblyContext& /*context*/ )
  {
    /*
    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.get_dof_indices(_flow_vars.u()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
    context.get_element_fe(_flow_vars.u())->get_JxW();

    const std::vector<std::vector<libMesh::Real> >& u_phi =
    context.get_element_fe(this->_flow_vars.u())->get_phi();

    // Get residuals
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(_flow_vars.u()); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(_flow_vars.v()); // R_{v}
    libMesh::DenseSubVector<libMesh::Number> *Fw = NULL;
    if(this->_flow_vars.dim() == 3)
    {
    Fw = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}
    }

    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.get_element_qrule().n_points();

    libMesh::FEBase* fe = context.get_element_fe(this->_flow_vars.u());

    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
    libMesh::RealGradient g = this->_flow_stab_helper.compute_g( fe, context, qp );
    libMesh::RealTensor G = this->_flow_stab_helper.compute_G( fe, context, qp );

    libMesh::RealGradient U( context.fixed_interior_value( this->_flow_vars.u(), qp ),
    context.fixed_interior_value( this->_flow_vars.v(), qp ) );
    if( this->_flow_vars.dim() == 3 )
    {
    U(2) = context.fixed_interior_value( this->_flow_vars.w(), qp );
    }

    libMesh::Real tau_E = this->_temp_stab_helper.compute_tau_energy( context, G, _rho, _Cp, _k,  U, false );

    libMesh::Real RE = this->_temp_stab_helper.compute_res_energy_transient( context, qp, _rho, _Cp );


    for (unsigned int i=0; i != n_u_dofs; i++)
    {
    Fu(i) += -_rho*_beta_T*tau_E*RE*_g(0)*u_phi[i][qp]*JxW[qp];

    Fv(i) += -_rho*_beta_T*tau_E*RE*_g(1)*u_phi[i][qp]*JxW[qp];

    if (this->_flow_vars.dim() == 3)
    {
    (*Fw)(i) += -_rho*_beta_T*tau_E*RE*_g(2)*u_phi[i][qp]*JxW[qp];
    }

    if (compute_jacobian)
    {
    libmesh_not_implemented();
    } // End compute_jacobian check

    } // End i dof loop
    } // End quadrature loop
    */
  }

  template<class Mu>
  void BoussinesqBuoyancySPGSMStabilization<Mu>::register_parameter
  ( const std::string & param_name,
    libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
    const
  {
    ParameterUser::register_parameter(param_name, param_pointer);
    _mu.register_parameter(param_name, param_pointer);
  }


} // namespace GRINS

// Instantiate
INSTANTIATE_INC_NS_SUBCLASS(BoussinesqBuoyancySPGSMStabilization);
