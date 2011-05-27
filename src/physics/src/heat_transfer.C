//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - a low Mach number Navier-Stokes Finite-Element Solver
//
// Copyright (C) 2010,2011 The PECOS Development Team
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
// Definitions for the HeatTransfer class.
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "heat_transfer.h"

void GRINS::HeatTransfer::read_input_options( GetPot& input )
{
  this->_FE_family =
    libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/HeatTrans/FE_family", "LAGRANGE") );

  this->_T_order =
    libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/HeatTransfer/T_order", "SECOND") );

  this->_rho = input("Physics/HeatTransfer/rho", 1.0); //TODO: same as Incompressible NS
  this->_Cp  = input("Physics/HeatTransfer/Cp", 1.0);
  this->_k  = input("Physics/HeatTransfer/k", 1.0);

  return;
}

void GRINS::HeatTransfer::init_variables( libMesh::FEMSystem* system )
{
  // Get libMesh to assign an index for each variable
  this->_dim = system->get_mesh().mesh_dimension();

  _T_var = system->add_variable( "T", this->_T_order, _FE_family);

  // Now build the local map
  this->build_local_variable_map();

  return;
}

void GRINS::HeatTransfer::build_local_variable_map()
{
  _var_map["T"] = _T_var;

  this->_local_variable_map_built = true;

  return;
}

void GRINS::HeatTransfer::set_time_evolving_vars( libMesh::FEMSystem* system )
{
  const unsigned int dim = system->get_mesh().mesh_dimension();

  // Tell the system to march temperature forward in time
  system->time_evolving(_T_var);

  return;
}

void GRINS::HeatTransfer::init_context( libMesh::DiffContext &context )
{
  libMesh::FEMContext &c = libmesh_cast_ref<libMesh::FEMContext&>(context);

  // We should prerequest all the data
  // we will need to build the linear system
  // or evaluate a quantity of interest.
  c.element_fe_var[_T_var]->get_JxW();
  c.element_fe_var[_T_var]->get_phi();
  c.element_fe_var[_T_var]->get_dphi();
  c.element_fe_var[_T_var]->get_xyz();

  c.side_fe_var[_T_var]->get_JxW();
  c.side_fe_var[_T_var]->get_phi();
  c.side_fe_var[_T_var]->get_dphi();
  c.side_fe_var[_T_var]->get_xyz();

  //TODO: _u_var is registered so can we assume things available for it in FEMContext?

  return;
}

bool GRINS::HeatTransfer::element_time_derivative( bool request_jacobian,
								 libMesh::DiffContext& context,
								 libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("HeatTransfer::element_time_derivative");
#endif

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // The number of local degrees of freedom in each variable.
  const unsigned int n_T_dofs = c.dof_indices_var[_T_var].size();

  //TODO: check n_T_dofs is same as n_u_dofs, n_v_dofs, n_w_dofs

  // We get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[_T_var]->get_JxW();

  // The temperature shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& T_phi =
    c.element_fe_var[_T_var]->get_phi();

  // The velocity shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& vel_phi =
    c.element_fe_var[_u_var]->get_phi();


  // The temperature shape function gradients (in global coords.)
  // at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& T_gblgradphivec =
    c.element_fe_var[_T_var]->get_dphi();

  // Physical location of the quadrature points
  const std::vector<libMesh::Point>& T_qpoint =
    c.element_fe_var[_T_var]->get_xyz();

  // The subvectors and submatrices we need to fill:
  //
  // K_{\alpha \beta} = R_{\alpha},{\beta} = \partial{ R_{\alpha} } / \partial{ {\beta} } (where R denotes residual)
  // e.g., for \alpha = T and \beta = v we get: K_{Tu} = R_{T},{u}
  //

  libMesh::DenseSubMatrix<Number> &KTT = *c.elem_subjacobians[_T_var][_T_var]; // R_{T},{T}

  libMesh::DenseSubMatrix<Number> &KTu = *c.elem_subjacobians[_T_var][_u_var]; // R_{T},{u}
  libMesh::DenseSubMatrix<Number> &KTv = *c.elem_subjacobians[_T_var][_v_var]; // R_{T},{v}
  libMesh::DenseSubMatrix<Number> &KTw = *c.elem_subjacobians[_T_var][_w_var]; // R_{T},{w}

  libMesh::DenseSubVector<Number> &FT = *c.elem_subresiduals[_T_var]; // R_{T}

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
      libMesh::Number T, u, v, w;
      T = c.interior_value(_T_var, qp);
      u = c.interior_value(_u_var, qp);
      v = c.interior_value(_v_var, qp);
      if (_dim == 3)
        w = c.interior_value(_w_var, qp);

      libMesh::Gradient gradTvec;
      gradTvec = c.interior_gradient(_T_var, qp);

      libMesh::NumberVectorValue Uvec (u,v);
      if (_dim == 3)
        Uvec(2) = w;

      // Value of the heat source function at this quadrature point.
      libMesh::Number heat_source = this->heat_source(T_qpoint[qp]);

      // First, an i-loop over the  degrees of freedom.
      for (unsigned int i=0; i != n_T_dofs; i++)
        {
          FT(i) += JxW[qp] *
                   (-_rho*_Cp*T_phi[i][qp]*(Uvec*gradTvec) // convection term
                    -_k*(T_gblgradphivec[i][qp]*gradTvec) // diffusion term
                    +T_phi[i][qp]*heat_source);           // source term

          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);

              for (unsigned int j=0; j != n_T_dofs; j++)
                {
                  // TODO: precompute some terms like:
                  // _rho*_Cp*T_phi[i][qp]*(vel_phi[j][qp]*T_gblgradphivec[j][qp])

                  KTT(i,j) += JxW[qp] *
                              (-_rho*_Cp*T_phi[i][qp]*(Uvec*T_gblgradphivec[j][qp])  // convection term
                               -_k*(T_gblgradphivec[i][qp]*T_gblgradphivec[j][qp])); // diffusion term
                } // end of the inner dof (j) loop

              // Matrix contributions for the Tu, Tv and Tw couplings (n_T_dofs same as n_u_dofs, n_v_dofs and n_w_dofs)
              for (unsigned int j=0; j != n_T_dofs; j++)
                {
                  KTu(i,j) += JxW[qp]*(-_rho*_Cp*T_phi[i][qp]*(vel_phi[j][qp]*T_gblgradphivec[j][qp](0)));
                  KTv(i,j) += JxW[qp]*(-_rho*_Cp*T_phi[i][qp]*(vel_phi[j][qp]*T_gblgradphivec[j][qp](1)));
                  if (_dim == 3)
                     KTw(i,j) += JxW[qp]*(-_rho*_Cp*T_phi[i][qp]*(vel_phi[j][qp]*T_gblgradphivec[j][qp](2)));
                } // end of the inner dof (j) loop

            } // end - if (request_jacobian && c.elem_solution_derivative)

        } // end of the outer dof (i) loop
    } // end of the quadrature point (qp) loop

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("HeatTransfer::element_time_derivative");
#endif

  return request_jacobian;
}

bool GRINS::HeatTransfer::element_constraint( bool request_jacobian,
							    libMesh::DiffContext& context,
							    libMesh::FEMSystem* system )
{
  return request_jacobian;
}


bool GRINS::HeatTransfer::side_time_derivative( bool request_jacobian,
							      libMesh::DiffContext& context,
							      libMesh::FEMSystem* system )
{
  return request_jacobian;
}

bool GRINS::HeatTransfer::side_constraint( bool request_jacobian,
							 libMesh::DiffContext& context,
							 libMesh::FEMSystem* system )
{
#ifdef USE_GRVY_TIMERS
  this->_timer->BeginTimer("HeatTransfer::side_constraint");
#endif

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // The number of local degrees of freedom in T variable.
  const unsigned int n_T_dofs = c.dof_indices_var[_T_var].size();

  // We get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weight for side integration.
  const std::vector<libMesh::Real> &JxW_side = c.side_fe_var[_T_var]->get_JxW();

  // The temperature shape functions at side quadrature points.
  const std::vector<std::vector<libMesh::Real> >& T_phi_side =
    c.side_fe_var[_T_var]->get_phi();

  // Physical location of the quadrature points on the side.
  const std::vector<libMesh::Point>& T_qpoint = c.side_fe_var[_T_var]->get_xyz();

  // The subvectors and submatrices we need to fill:
  //
  // KTT, FT
  //

  libMesh::DenseSubMatrix<Number> &KTT = *c.elem_subjacobians[_T_var][_T_var]; // R_{T},{T}

  libMesh::DenseSubVector<Number> &FT = *c.elem_subresiduals[_T_var]; // R_{T}

  // For this example we will use Dirichlet temperature boundary
  // conditions imposed at each timestep via the penalty method.

  // The penalty value.  \f$ \frac{1}{\epsilon} \f$
  const libMesh::Real T_penalty = 1.e10;

  const short int boundary_id =
    system->get_mesh().boundary_info->boundary_id(c.elem, c.side);
  libmesh_assert (boundary_id != libMesh::BoundaryInfo::invalid_id);

  unsigned int n_sidepoints = c.side_qrule->n_points();
  for (unsigned int qp=0; qp != n_sidepoints; qp++)
    {
      // Compute the solution at the old Newton iterate
      libMesh::Number T = c.side_value(_T_var, qp);

      // TODO: make it more general
      short int left_id, right_id, top_id, bottom_id, front_id, back_id;
      if (_dim == 1)
        {
          left_id  = 0; // x=0
          right_id = 1; // x=1
        }
      if (_dim == 2)
        {
          left_id   = 3; // x=0
          right_id  = 1; // x=1
          bottom_id = 0; // y=0
          top_id    = 2; // y=1
        }
      if (_dim == 3)
        {
          left_id   = 4; // x=0
          right_id  = 2; // x=1
          bottom_id = 1; // y=0
          top_id    = 3; // y=1
          back_id   = 0; // z=0
          front_id  = 5; // z=1
        }

      if ((boundary_id == bottom_id) || (boundary_id == top_id))
        break;

      // Boundary data coming from true solution
      libMesh::Real T_value, T_value_left = 0.0, T_value_right = 1.0;

      // set T=0 on the left.
      if ((boundary_id == left_id))
        T_value = T_value_left;
      // set T=1 on the right.
      if ((boundary_id == right_id))
        T_value = T_value_right;

      for (unsigned int i=0; i != n_T_dofs; i++)
        {
          FT(i) += JxW_side[qp] * T_penalty *
                   (T - T_value) * T_phi_side[i][qp];

          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);

              for (unsigned int j=0; j != n_T_dofs; j++)
                {
                  KTT(i,j) += JxW_side[qp] * T_penalty *
                              T_phi_side[i][qp] * T_phi_side[j][qp];
                } // end of the inner dof (j) loop

            } // end - if (request_jacobian && c.elem_solution_derivative)

        } // end of the outer dof (i) loop
    } // end of the quadrature point (qp) loop

#ifdef USE_GRVY_TIMERS
  this->_timer->EndTimer("HeatTransfer::side_constraint");
#endif

  return request_jacobian;
}

bool GRINS::HeatTransfer::mass_residual( bool request_jacobian,
						       libMesh::DiffContext& context,
						       libMesh::FEMSystem* system )
{
  // TODO: account for 'rho*Cp' factor in mass matrix
  return request_jacobian;
}

libMesh::Number GRINS::HeatTransfer::heat_source(const libMesh::Point &pt_xyz)
{
  // TODO: add heat source options
  return libMesh::Number(0.);
}
