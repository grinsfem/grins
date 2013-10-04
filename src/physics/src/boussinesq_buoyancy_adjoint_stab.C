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


// This class
#include "grins_config.h"
#include "grins/boussinesq_buoyancy_adjoint_stab.h"

// GRINS
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace GRINS
{

  BoussinesqBuoyancyAdjointStabilization::BoussinesqBuoyancyAdjointStabilization( const std::string& physics_name, const GetPot& input )
    : BoussinesqBuoyancyBase(physics_name,input),
      /* \todo Do we want to have these come from a BoussinesqBuoyancyAdjointStabilization section instead? */
      _rho( input("Physics/"+incompressible_navier_stokes+"/rho", 1.0) ),
      _mu( input("Physics/"+incompressible_navier_stokes+"/mu", 1.0) ),
      _stab_helper( input ),
      _p_var_name( input("Physics/VariableNames/pressure", p_var_name_default ) ),
      _P_FE_family( libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+incompressible_navier_stokes+"/FE_family", "LAGRANGE") ) ),
      _P_order( libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+incompressible_navier_stokes+"/P_order", "FIRST") ) )
  {
    return;
  }

  BoussinesqBuoyancyAdjointStabilization::~BoussinesqBuoyancyAdjointStabilization()
  {
    return;
  }

  void BoussinesqBuoyancyAdjointStabilization::init_context( AssemblyContext& context )
  {
    context.get_element_fe(this->_p_var)->get_dphi();

    context.get_element_fe(this->_u_var)->get_dphi();
    context.get_element_fe(this->_u_var)->get_d2phi();

    return;
  }

  void BoussinesqBuoyancyAdjointStabilization::init_variables( libMesh::FEMSystem* system )
  {
    // First call base class
    BoussinesqBuoyancyBase::init_variables(system);

    _p_var = system->add_variable( _p_var_name, this->_P_order, _P_FE_family);

    return;
  }

  void BoussinesqBuoyancyAdjointStabilization::element_time_derivative( bool compute_jacobian,
                                                                        AssemblyContext& context,
                                                                        CachedValues& /*cache*/ )
  {
#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->BeginTimer("BoussinesqBuoyancyAdjointStabilization::element_time_derivative");
#endif

    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.get_dof_indices(_u_var).size();
    const unsigned int n_p_dofs = context.get_dof_indices(_p_var).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(_u_var)->get_JxW();

    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_u_var)->get_dphi();

    const std::vector<std::vector<libMesh::RealTensor> >& u_hessphi =
      context.get_element_fe(this->_u_var)->get_d2phi();

    const std::vector<std::vector<libMesh::RealGradient> >& p_dphi =
      context.get_element_fe(this->_p_var)->get_dphi();

    // Get residuals
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(_u_var); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(_v_var); // R_{v}
    libMesh::DenseSubVector<libMesh::Number> *Fw = NULL;
    if(this->_dim == 3)
      Fw = &context.get_elem_residual(this->_w_var); // R_{w}

    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(this->_p_var); // R_{p}

    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.get_element_qrule().n_points();

    libMesh::FEBase* fe = context.get_element_fe(this->_u_var);

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::RealGradient U( context.interior_value( this->_u_var, qp ),
                                 context.interior_value( this->_v_var, qp ) );
        if( this->_dim == 3 )
          {
            U(2) = context.interior_value( this->_w_var, qp );
          }

        libMesh::Real tau_M = this->_stab_helper.compute_tau_momentum( context, qp, g, G, this->_rho, U, this->_mu, this->_is_steady );

        // Compute the solution & its gradient at the old Newton iterate.
        libMesh::Number T;
        T = context.interior_value(_T_var, qp);

        libMesh::RealGradient residual = -_rho_ref*_beta_T*(T-_T_ref)*_g;

        // First, an i-loop over the velocity degrees of freedom.
        // We know that n_u_dofs == n_v_dofs so we can compute contributions
        // for both at the same time.
        for (unsigned int i=0; i != n_p_dofs; i++)
          {
            Fp(i) += tau_M*residual*p_dphi[i][qp]*JxW[qp];
          }

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) += ( this->_rho*U*u_gradphi[i][qp]
                       + this->_mu*( u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) + u_hessphi[i][qp](2,2) ) )*tau_M*residual(0)*JxW[qp];

            Fv(i) += ( this->_rho*U*u_gradphi[i][qp]
                       + this->_mu*( u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) + u_hessphi[i][qp](2,2) ) )*tau_M*residual(1)*JxW[qp];

            if (_dim == 3)
              {
                (*Fw)(i) += ( this->_rho*U*u_gradphi[i][qp]
                              + this->_mu*( u_hessphi[i][qp](0,0) + u_hessphi[i][qp](1,1) + u_hessphi[i][qp](2,2) ) )*tau_M*residual(2)*JxW[qp];
              }

            if (compute_jacobian)
              {
                libmesh_not_implemented();
              } // End compute_jacobian check

          } // End i dof loop
      } // End quadrature loop

#ifdef GRINS_USE_GRVY_TIMERS
    this->_timer->EndTimer("BoussinesqBuoyancyAdjointStabilization::element_time_derivative");
#endif

    return;
  }

} // namespace GRINS
