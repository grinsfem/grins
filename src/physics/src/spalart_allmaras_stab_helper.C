//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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
#include "grins/spalart_allmaras_stab_helper.h"

//libMesh
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  SpalartAllmarasStabilizationHelper::SpalartAllmarasStabilizationHelper(const GetPot& input)
    : StabilizationHelper(),
      _C( input("Stabilization/tau_constant_vel", input("Stabilization/tau_constant", 1 ) ) ),
      _tau_factor( input("Stabilization/tau_factor_vel", input("Stabilization/tau_factor", 0.5 ) ) ),
      _flow_vars(input),
      _turbulence_vars(input),
      _spalart_allmaras_helper(input)
  {
    return;
  }

  SpalartAllmarasStabilizationHelper::~SpalartAllmarasStabilizationHelper()
  {
    return;
  }

  void SpalartAllmarasStabilizationHelper::init( libMesh::FEMSystem& system )
  {
    this->_flow_vars.init(&system);
    this->_turbulence_vars.init(&system);

    // Init the variables belonging to SA helper
    _spalart_allmaras_helper.init_variables(&system);

    this->_dim = system.get_mesh().mesh_dimension();

    return;
  }

  libMesh::Real SpalartAllmarasStabilizationHelper::compute_res_spalart_steady( AssemblyContext& context,
                                                                                unsigned int qp, const libMesh::Real rho, const libMesh::Real mu, const libMesh::Real distance_qp ) const
  {
    // The flow velocity
    libMesh::Number u,v;
    u = context.interior_value(this->_flow_vars.u_var(), qp);
    v = context.interior_value(this->_flow_vars.v_var(), qp);

    libMesh::NumberVectorValue U(u,v);
    if ( context.get_system().get_mesh().mesh_dimension() == 3 )
      U(2) = context.interior_value(this->_flow_vars.w_var(), qp);

    libMesh::RealGradient grad_u = context.fixed_interior_gradient(this->_flow_vars.u_var(), qp);
    libMesh::RealGradient grad_v = context.fixed_interior_gradient(this->_flow_vars.v_var(), qp);

    libMesh::Number nu_value = context.interior_value(this->_turbulence_vars.nu_var(), qp);

    libMesh::RealGradient grad_nu = context.fixed_interior_gradient(this->_turbulence_vars.nu_var(), qp);

    libMesh::RealTensor hess_nu = context.fixed_interior_hessian(this->_turbulence_vars.nu_var(), qp);

    // The convection term
    libMesh::Number rhoUdotGradnu = rho*(U*grad_nu);

    // The diffusion term
    libMesh::Number inv_sigmadivnuplusnuphysicalGradnu = (1./this->_spalart_allmaras_helper._sigma)*(grad_nu*grad_nu + ((nu_value + mu)*(hess_nu(0,0) + hess_nu(1,1) + (this->_dim == 3)?hess_nu(2,2):0)) + this->_spalart_allmaras_helper._cb2*grad_nu*grad_nu);

    // The source term
    libMesh::Real _vorticity_value_qp = this->_spalart_allmaras_helper._vorticity(context, qp);
    libMesh::Real _S_tilde = this->_spalart_allmaras_helper._source_fn(nu_value, mu, distance_qp, _vorticity_value_qp);
    libMesh::Real source_term = this->_spalart_allmaras_helper._cb1*_S_tilde*nu_value;

    // The destruction term
    libMesh::Real _fw = this->_spalart_allmaras_helper._destruction_fn(nu_value, distance_qp, _S_tilde);
    libMesh::Real destruction_term =  this->_spalart_allmaras_helper._cw1*_fw*pow(nu_value/distance_qp, 2.);

    return rhoUdotGradnu + source_term + inv_sigmadivnuplusnuphysicalGradnu - destruction_term;
  }

  void SpalartAllmarasStabilizationHelper::compute_res_spalart_steady_and_derivs
  ( AssemblyContext& /*context*/,
    unsigned int /*qp*/, const libMesh::Real /*rho*/, const libMesh::Real /*mu*/,
    libMesh::Gradient& /*res_M*/,
    libMesh::Tensor&   /*d_res_M_dgradp*/,
    libMesh::Tensor&   /*d_res_M_dU*/,
    libMesh::Gradient& /*d_res_Muvw_dgraduvw*/,
    libMesh::Tensor&   /*d_res_Muvw_dhessuvw*/
    ) const
  {
    // To be filled when we start using analytic jacobians with SA
    libmesh_not_implemented();
  }


  libMesh::Real SpalartAllmarasStabilizationHelper::compute_res_spalart_transient( AssemblyContext& context, unsigned int qp, const libMesh::Real rho ) const
  {
    libMesh::Number nu_dot = context.interior_value(this->_turbulence_vars.nu_var(), qp);

    return rho*nu_dot;
  }


  void SpalartAllmarasStabilizationHelper::compute_res_spalart_transient_and_derivs
  ( AssemblyContext& /*context*/,
    unsigned int /*qp*/,
    const libMesh::Real /*rho*/,
    libMesh::RealGradient& /*res_M*/,
    libMesh::Real& /*d_res_Muvw_duvw*/
    ) const
  {
    libmesh_not_implemented();
  }

} // namespace GRINS
