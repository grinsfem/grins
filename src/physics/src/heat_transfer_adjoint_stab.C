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

#include "heat_transfer_adjoint_stab.h"

GRINS::HeatTransferAdjointStabilization::HeatTransferAdjointStabilization( const std::string& physics_name, 
									   const GetPot& input )
  : GRINS::HeatTransferStabilizationBase(physics_name,input)
{
  this->read_input_options(input);

  return;
}

GRINS::HeatTransferAdjointStabilization::~HeatTransferAdjointStabilization()
{
  return;
}

void GRINS::HeatTransferAdjointStabilization::read_input_options( const GetPot& )
{
  return;
}

bool GRINS::HeatTransferAdjointStabilization::element_time_derivative( bool request_jacobian,
								       libMesh::DiffContext& context,
								       libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("HeatTransferAdjointStabilization::element_time_derivative");
#endif

  libMesh::FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // The number of local degrees of freedom in each variable.
  const unsigned int n_T_dofs = c.dof_indices_var[this->_T_var].size();

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[this->_T_var]->get_JxW();

  const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
    c.element_fe_var[this->_T_var]->get_dphi();

  const std::vector<std::vector<libMesh::RealTensor> >& T_hessphi =
    c.element_fe_var[this->_T_var]->get_d2phi();

  libMesh::DenseSubVector<Number> &FT = *c.elem_subresiduals[this->_T_var]; // R_{T}

  unsigned int n_qpoints = c.element_qrule->n_points();

  bool is_steady = (system->time_solver)->is_steady();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      libMesh::FEBase* fe = c.element_fe_var[this->_T_var];

      libMesh::RealGradient g = this->_stab_helper.compute_g( fe, c, qp );
      libMesh::RealTensor G = this->_stab_helper.compute_G( fe, c, qp );

      libMesh::RealGradient U( c.interior_value( this->_u_var, qp ),
			       c.interior_value( this->_v_var, qp ) );
      if( this->_dim == 3 )
	U(2) = c.interior_value( this->_w_var, qp );
      
      libMesh::Real tau_E = this->_stab_helper.compute_tau_energy( c, G, _rho, _Cp, _k,  U, is_steady );

      libMesh::Real RE_s = this->compute_res_steady( c, qp );

      for (unsigned int i=0; i != n_T_dofs; i++)
        {
          FT(i) += tau_E*RE_s*( _rho*_Cp*U*T_gradphi[i][qp]
				+ _k*(T_hessphi[i][qp](0,0) + T_hessphi[i][qp](1,1) + T_hessphi[i][qp](2,2)) 
				)*JxW[qp];
	}

    }

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("HeatTransferAdjointStabilization::element_time_derivative");
#endif
  return request_jacobian;
}

bool GRINS::HeatTransferAdjointStabilization::mass_residual( bool request_jacobian,
							     libMesh::DiffContext& context,
							     libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("HeatTransferAdjointStabilization::mass_residual");
#endif

  libMesh::FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // The number of local degrees of freedom in each variable.
  const unsigned int n_T_dofs = c.dof_indices_var[this->_T_var].size();

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[this->_T_var]->get_JxW();

  const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
    c.element_fe_var[this->_T_var]->get_dphi();

  const std::vector<std::vector<libMesh::RealTensor> >& T_hessphi =
    c.element_fe_var[this->_T_var]->get_d2phi();

  libMesh::DenseSubVector<Number> &FT = *c.elem_subresiduals[this->_T_var]; // R_{T}

  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      libMesh::FEBase* fe = c.element_fe_var[this->_T_var];

      libMesh::RealGradient g = this->_stab_helper.compute_g( fe, c, qp );
      libMesh::RealTensor G = this->_stab_helper.compute_G( fe, c, qp );

      libMesh::RealGradient U( c.fixed_interior_value( this->_u_var, qp ),
			       c.fixed_interior_value( this->_v_var, qp ) );
      if( this->_dim == 3 )
	U(2) = c.fixed_interior_value( this->_w_var, qp );
      
      libMesh::Real tau_E = this->_stab_helper.compute_tau_energy( c, G, _rho, _Cp, _k,  U, false );

      libMesh::Real RE_t = this->compute_res_transient( c, qp );

      for (unsigned int i=0; i != n_T_dofs; i++)
        {
          FT(i) -= tau_E*RE_t*( _rho*_Cp*U*T_gradphi[i][qp]
				+ _k*(T_hessphi[i][qp](0,0) + T_hessphi[i][qp](1,1) + T_hessphi[i][qp](2,2)) 
				)*JxW[qp];
	}

    }

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("HeatTransferAdjointStabilization::mass_residual");
#endif
  return request_jacobian;
}

bool GRINS::HeatTransferAdjointStabilization::element_constraint( bool request_jacobian,
										libMesh::DiffContext& context,
										libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  //this->_timer->BeginTimer("HeatTransferAdjointStabilization::element_constraint");
#endif

  //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

#ifdef USE_GRVY_TIMERS
  //this->_timer->EndTimer("HeatTransferAdjointStabilization::element_constraint");
#endif

  return request_jacobian;
}

bool GRINS::HeatTransferAdjointStabilization::side_time_derivative( bool request_jacobian,
										  libMesh::DiffContext& context,
										  libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
      //this->_timer->BeginTimer("HeatTransferAdjointStabilization::side_time_derivative");
#endif
      //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

#ifdef USE_GRVY_TIMERS
      //this->_timer->EndTimer("HeatTransferAdjointStabilization::side_time_derivative");
#endif

  return request_jacobian;
}

bool GRINS::HeatTransferAdjointStabilization::side_constraint( bool request_jacobian,
									     libMesh::DiffContext& context,
									     libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  //this->_timer->BeginTimer("HeatTransferAdjointStabilization::side_constraint");
#endif

  //FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

#ifdef USE_GRVY_TIMERS
  //this->_timer->EndTimer("HeatTransferAdjointStabilization::side_constraint");
#endif

  return request_jacobian;
}
