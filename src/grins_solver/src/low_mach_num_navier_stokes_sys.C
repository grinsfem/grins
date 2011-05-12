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
// Definitions for the LowMachNumberNavierStokesSystem class.
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "low_mach_num_navier_stokes_sys.h"

void GRINS::LowMachNumberNavierStokesSystem::init_data()
{
  const unsigned int dim = this->get_mesh().mesh_dimension();

  const libMeshEnums::FEFamily FE_family = libMeshEnums::LAGRANGE;
  const libMeshEnums::Order V_order = libMeshEnums::SECOND;
  const libMeshEnums::Order P_order = static_cast<libMeshEnums::Order>(V_order-1);

  // Flag for plain/standard incompressible flow
  const bool plain_ic_flag = true;

  u_var = this->add_variable( "u", V_order, FE_family);
  v_var = this->add_variable( "v", V_order, FE_family);
  if (dim == 3)
    w_var = this->add_variable( "w", V_order, FE_family);
  p_var = this->add_variable( "p", P_order, FE_family);

  libmesh_assert (plain_ic_flag == true);

  // Do the parent's initialization after variables are defined
  FEMSystem::init_data();

  // Tell the system to march velocity forward in time, but
  // leave p as a constraint only
  this->time_evolving(u_var);
  this->time_evolving(v_var);
  if (dim == 3)
    this->time_evolving(w_var);

  return;
}

void GRINS::LowMachNumberNavierStokesSystem::init_context(
						libMesh::DiffContext &context )
{
  libMesh::FEMContext &c = libmesh_cast_ref<libMesh::FEMContext&>(context);

  // We should prerequest all the data
  // we will need to build the linear system
  // or evaluate a quantity of interest.
  c.element_fe_var[u_var]->get_JxW();
  c.element_fe_var[u_var]->get_phi();
  c.element_fe_var[u_var]->get_dphi();
  c.element_fe_var[u_var]->get_xyz();

  c.element_fe_var[p_var]->get_phi();
  c.element_fe_var[p_var]->get_xyz();

  c.side_fe_var[u_var]->get_JxW();
  c.side_fe_var[u_var]->get_phi();
  c.side_fe_var[u_var]->get_dphi();
  c.side_fe_var[u_var]->get_xyz();

  return;
}

bool GRINS::LowMachNumberNavierStokesSystem::element_time_derivative(
						bool request_jacobian,
						libMesh::DiffContext& context )
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  const unsigned int dim = this->get_mesh().mesh_dimension();

  // The number of local degrees of freedom in each variable.
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
  const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();

  // Check number of dofs is same for u_var, v_var and w_var.
  libmesh_assert (n_u_dofs == c.dof_indices_var[v_var].size());
  if (dim == 3)
    libmesh_assert (n_u_dofs == c.dof_indices_var[w_var].size());

  // We get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[u_var]->get_JxW();

  // The velocity shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& vel_phi =
    c.element_fe_var[u_var]->get_phi();

  // The velocity shape function gradients (in global coords.)
  // at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& vel_gblgradphivec =
    c.element_fe_var[u_var]->get_dphi();

  // The pressure shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& p_phi =
    c.element_fe_var[p_var]->get_phi();

  // Physical location of the quadrature points
  const std::vector<libMesh::Point>& u_qpoint =
    c.element_fe_var[u_var]->get_xyz();

  // The subvectors and submatrices we need to fill:
  //
  // K_{\alpha \beta} = R_{\alpha},{\beta} = \partial{ R_{\alpha} } / \partial{ {\beta} } )where R is residual)
  // e.g., for \alpha = v and \beta = u we get: Kvu = Rv,u
  // Note that Kpu, Kpv, Kpw and Fp comes as constraint.
  //
  if (dim != 3)
    w_var = u_var; // for convenience

  libMesh::DenseSubMatrix<Number> &Kuu = *c.elem_subjacobians[u_var][u_var]; // Ru,u
  libMesh::DenseSubMatrix<Number> &Kuv = *c.elem_subjacobians[u_var][v_var]; // Ru,v
  libMesh::DenseSubMatrix<Number> &Kuw = *c.elem_subjacobians[u_var][w_var]; // Ru,w

  libMesh::DenseSubMatrix<Number> &Kvu = *c.elem_subjacobians[v_var][u_var]; // Rv,u
  libMesh::DenseSubMatrix<Number> &Kvv = *c.elem_subjacobians[v_var][v_var]; // Rv,v
  libMesh::DenseSubMatrix<Number> &Kvw = *c.elem_subjacobians[v_var][w_var]; // Rv,w

  libMesh::DenseSubMatrix<Number> &Kwu = *c.elem_subjacobians[w_var][u_var]; // Rw,u
  libMesh::DenseSubMatrix<Number> &Kwv = *c.elem_subjacobians[w_var][v_var]; // Rw,v
  libMesh::DenseSubMatrix<Number> &Kww = *c.elem_subjacobians[w_var][w_var]; // Rw,w

  libMesh::DenseSubMatrix<Number> &Kup = *c.elem_subjacobians[u_var][p_var]; // Ru,p
  libMesh::DenseSubMatrix<Number> &Kvp = *c.elem_subjacobians[v_var][p_var]; // Ru,p
  libMesh::DenseSubMatrix<Number> &Kwp = *c.elem_subjacobians[w_var][p_var]; // Rw,p

  libMesh::DenseSubVector<Number> &Fu = *c.elem_subresiduals[u_var]; // Ru
  libMesh::DenseSubVector<Number> &Fv = *c.elem_subresiduals[v_var]; // Rv
  libMesh::DenseSubVector<Number> &Fw = *c.elem_subresiduals[w_var]; // Rw

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
      p = c.interior_value(p_var, qp);
      u = c.interior_value(u_var, qp);
      v = c.interior_value(v_var, qp);
      if (dim == 3)
        w = c.interior_value(w_var, qp);

      libMesh::Gradient graduvec, gradvvec, gradwvec;
      graduvec = c.interior_gradient(u_var, qp);
      gradvvec = c.interior_gradient(v_var, qp);
      if (dim == 3)
        gradwvec = c.interior_gradient(w_var, qp);

      libMesh::NumberVectorValue Uvec (u,v);
      if (dim == 3)
        Uvec(2) = w;

      const libMesh::Number  graduvec_x = graduvec(0);
      const libMesh::Number  graduvec_y = graduvec(1);
      const libMesh::Number  graduvec_z = (dim == 3)?graduvec(2):0;
      const libMesh::Number  gradvvec_x = gradvvec(0);
      const libMesh::Number  gradvvec_y = gradvvec(1);
      const libMesh::Number  gradvvec_z = (dim == 3)?gradvvec(2):0;
      const libMesh::Number  gradwvec_x = (dim == 3)?gradwvec(0):0;
      const libMesh::Number  gradwvec_y = (dim == 3)?gradwvec(1):0;
      const libMesh::Number  gradwvec_z = (dim == 3)?gradwvec(2):0;

      // Value of the forcing function at this quadrature point.
      libMesh::Point force = this->forcing(u_qpoint[qp]);

      const libMesh::Number rho = 1.0;           // TODO: set these properly!!!
      const libMesh::Number mu  = 1.0e-0;        // TODO: set these properly!!!

      // First, an i-loop over the velocity degrees of freedom.
      // We know that n_u_dofs == n_v_dofs so we can compute contributions
      // for both at the same time.
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fu(i) += JxW[qp] *
                   (-rho*vel_phi[i][qp]*(Uvec*graduvec)     // convection term
                    +p*vel_gblgradphivec[i][qp](0)          // pressure term
                    -mu*(vel_gblgradphivec[i][qp]*graduvec) // diffusion term
                    +vel_phi[i][qp]*force(0));              // forcing function

          Fv(i) += JxW[qp] *
                   (-rho*vel_phi[i][qp]*(Uvec*gradvvec)     // convection term
                    +p*vel_gblgradphivec[i][qp](1)          // pressure term
                    -mu*(vel_gblgradphivec[i][qp]*gradvvec) // diffusion term
                    +vel_phi[i][qp]*force(1));              // forcing function
          if (dim == 3)
            {
              Fw(i) += JxW[qp] *
                       (-rho*vel_phi[i][qp]*(Uvec*gradwvec)     // convection term
                        +p*vel_gblgradphivec[i][qp](2)          // pressure term
                        -mu*(vel_gblgradphivec[i][qp]*gradwvec) // diffusion term
                        +vel_phi[i][qp]*force(2));              // forcing function
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
                              (-rho*vel_phi[i][qp]*(Uvec*vel_gblgradphivec[j][qp])       // convection term
                               -rho*vel_phi[i][qp]*graduvec_x*vel_phi[j][qp]             // convection term
                               -mu*(vel_gblgradphivec[i][qp]*vel_gblgradphivec[j][qp])); // diffusion term
                  Kuv(i,j) += JxW[qp] *
                              (-rho*vel_phi[i][qp]*graduvec_y*vel_phi[j][qp]);           // convection term

                  Kvv(i,j) += JxW[qp] *
                              (-rho*vel_phi[i][qp]*(Uvec*vel_gblgradphivec[j][qp])       // convection term
                               -rho*vel_phi[i][qp]*gradvvec_y*vel_phi[j][qp]             // convection term
                               -mu*(vel_gblgradphivec[i][qp]*vel_gblgradphivec[j][qp])); // diffusion term
                  Kvu(i,j) += JxW[qp] *
                              (-rho*vel_phi[i][qp]*gradvvec_x*vel_phi[j][qp]);           // convection term

                  if (dim == 3)
                    {
                      Kuw(i,j) += JxW[qp] *
                                  (-rho*vel_phi[i][qp]*graduvec_z*vel_phi[j][qp]);           // convection term

                      Kvw(i,j) += JxW[qp] *
                                  (-rho*vel_phi[i][qp]*gradvvec_z*vel_phi[j][qp]);           // convection term

                      Kww(i,j) += JxW[qp] *
                                  (-rho*vel_phi[i][qp]*(Uvec*vel_gblgradphivec[j][qp])       // convection term
                                   -rho*vel_phi[i][qp]*gradwvec_z*vel_phi[j][qp]             // convection term
                                   -mu*(vel_gblgradphivec[i][qp]*vel_gblgradphivec[j][qp])); // diffusion term
                      Kwu(i,j) += JxW[qp] *
                                  (-rho*vel_phi[i][qp]*gradwvec_x*vel_phi[j][qp]);           // convection term
                      Kwv(i,j) += JxW[qp] *
                                  (-rho*vel_phi[i][qp]*gradwvec_y*vel_phi[j][qp]);           // convection term
                    }
                } // end of the inner dof (j) loop

              // Matrix contributions for the up and vp couplings
              for (unsigned int j=0; j != n_p_dofs; j++)
                {
                  Kup(i,j) += JxW[qp]*vel_gblgradphivec[i][qp](0)*p_phi[j][qp];
                  Kvp(i,j) += JxW[qp]*vel_gblgradphivec[i][qp](1)*p_phi[j][qp];
                  if (dim == 3)
                    Kwp(i,j) += JxW[qp]*vel_gblgradphivec[i][qp](2)*p_phi[j][qp];
                } // end of the inner dof (j) loop

            } // end - if (request_jacobian && c.elem_solution_derivative)

        } // end of the outer dof (i) loop
    } // end of the quadrature point (qp) loop

  return request_jacobian;
}

bool GRINS::LowMachNumberNavierStokesSystem::element_constraint(
						bool request_jacobian,
						libMesh::DiffContext& context )
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  const unsigned int dim = this->get_mesh().mesh_dimension();

  // The number of local degrees of freedom in each variable.
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
  const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();

  // We get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration.
  const std::vector<libMesh::Real> &JxW =
    c.element_fe_var[u_var]->get_JxW();

  // The velocity shape function gradients (in global coords.)
  // at interior quadrature points.
  const std::vector<std::vector<libMesh::RealGradient> >& vel_gblgradphivec =
    c.element_fe_var[u_var]->get_dphi();

  // The pressure shape functions at interior quadrature points.
  const std::vector<std::vector<libMesh::Real> >& p_phi =
    c.element_fe_var[p_var]->get_phi();

  // The subvectors and submatrices we need to fill:
  //
  // Kpu, Kpv, Kpw, Fp
  //
  if (dim != 3)
    w_var = u_var; // for convenience

  libMesh::DenseSubMatrix<Number> &Kpu = *c.elem_subjacobians[p_var][u_var]; // Rp,u
  libMesh::DenseSubMatrix<Number> &Kpv = *c.elem_subjacobians[p_var][v_var]; // Rp,u
  libMesh::DenseSubMatrix<Number> &Kpw = *c.elem_subjacobians[p_var][w_var]; // Rp,u

  libMesh::DenseSubVector<Number> &Fp = *c.elem_subresiduals[p_var]; // Rp

  // Add the constraint given by the continuity equation.
  unsigned int n_qpoints = c.element_qrule->n_points();
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the velocity gradient at the old Newton iterate.
      libMesh::Gradient graduvec, gradvvec, gradwvec;
      graduvec = c.interior_gradient(u_var, qp);
      gradvvec = c.interior_gradient(v_var, qp);
      if (dim == 3)
        gradwvec = c.interior_gradient(w_var, qp);

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) += JxW[qp] * p_phi[i][qp] *
                   (graduvec(0) + gradvvec(1));
          if (dim == 3)
            Fp(i) += JxW[qp] * p_phi[i][qp] *
                     (gradwvec(2));

          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);

              for (unsigned int j=0; j != n_u_dofs; j++)
                {
                  Kpu(i,j) += JxW[qp]*p_phi[i][qp]*vel_gblgradphivec[j][qp](0);
                  Kpv(i,j) += JxW[qp]*p_phi[i][qp]*vel_gblgradphivec[j][qp](1);
                  if (dim == 3)
                    Kpw(i,j) += JxW[qp]*p_phi[i][qp]*vel_gblgradphivec[j][qp](2);
                } // end of the inner dof (j) loop

            } // end - if (request_jacobian && c.elem_solution_derivative)

        } // end of the outer dof (i) loop
    } // end of the quadrature point (qp) loop

  return request_jacobian;
}


bool GRINS::LowMachNumberNavierStokesSystem::side_time_derivative(
						bool request_jacobian,
						libMesh::DiffContext& context )
{
  return request_jacobian;
}

bool GRINS::LowMachNumberNavierStokesSystem::side_constraint(
						bool request_jacobian,
						libMesh::DiffContext& context )
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  const unsigned int dim = this->get_mesh().mesh_dimension();

  // The number of local degrees of freedom in u variable.
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();

  // We get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weight for side integration.
  const std::vector<libMesh::Real> &JxW_side = c.side_fe_var[u_var]->get_JxW();

  // The velocity shape functions at side quadrature points.
  const std::vector<std::vector<libMesh::Real> >& vel_phi_side =
    c.side_fe_var[u_var]->get_phi();

  // Physical location of the quadrature points on the side.
  const std::vector<libMesh::Point>& u_qpoint = c.side_fe_var[u_var]->get_xyz();

  // The subvectors and submatrices we need to fill:
  //
  // Kuu, Kvv, Kww, Fu, Fv, Fw
  //
  if (dim != 3)
    w_var = u_var; // for convenience

  libMesh::DenseSubMatrix<Number> &Kuu = *c.elem_subjacobians[u_var][u_var]; // Ru,u
  libMesh::DenseSubMatrix<Number> &Kvv = *c.elem_subjacobians[v_var][v_var]; // Rv,v
  libMesh::DenseSubMatrix<Number> &Kww = *c.elem_subjacobians[w_var][w_var]; // Rw,w

  libMesh::DenseSubVector<Number> &Fu = *c.elem_subresiduals[u_var]; // Ru
  libMesh::DenseSubVector<Number> &Fv = *c.elem_subresiduals[v_var]; // Rv
  libMesh::DenseSubVector<Number> &Fw = *c.elem_subresiduals[w_var]; // Rw

  // For this example we will use Dirichlet velocity boundary
  // conditions imposed at each timestep via the penalty method.

  // The penalty value.  \f$ \frac{1}{\epsilon} \f$
  const libMesh::Real vel_penalty = 1.e10;

  const short int boundary_id =
    this->get_mesh().boundary_info->boundary_id(c.elem, c.side);
  libmesh_assert (boundary_id != libMesh::BoundaryInfo::invalid_id);

  unsigned int n_sidepoints = c.side_qrule->n_points();
  for (unsigned int qp=0; qp != n_sidepoints; qp++)
    {
      // Compute the solution at the old Newton iterate
      libMesh::Number u = c.side_value(u_var, qp),
                      v = c.side_value(v_var, qp),
                      w = c.side_value(w_var, qp);

      // TODO: make it more general
      const short int left_id = (dim==3) ? 5 : 3; // TODO: CHECK THIS!!!
      const short int right_id = (dim==3) ? 5 : 1; // TODO: CHECK THIS!!!
      const short int top_id = (dim==3) ? 5 : 2; // TODO: CHECK THIS!!!

      // if (boundary_id == right_id)
      //  break; // its outflow

      // Boundary data coming from true solution
      libMesh::Real u_value = 0., v_value = 0., w_value = 0.;

      // For lid-driven cavity, set u=1 on the lid.
      if ((boundary_id == top_id))
        u_value = 1.;

      /* if (boundary_id == left_id)
        {
          const libMesh::Real x = u_qpoint[qp](0);
          const libMesh::Real y = u_qpoint[qp](1);
          const libMesh::Real z = u_qpoint[qp](2);
          u_value = 4.*y*(1.-y);
        } */

      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          Fu(i) += JxW_side[qp] * vel_penalty *
                   (u - u_value) * vel_phi_side[i][qp];
          Fv(i) += JxW_side[qp] * vel_penalty *
                   (v - v_value) * vel_phi_side[i][qp];
          if (dim == 3)
            Fw(i) += JxW_side[qp] * vel_penalty *
                     (w - w_value) * vel_phi_side[i][qp];

          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);

              for (unsigned int j=0; j != n_u_dofs; j++)
                {
                  Kuu(i,j) += JxW_side[qp] * vel_penalty *
                              vel_phi_side[i][qp] * vel_phi_side[j][qp];
                  Kvv(i,j) += JxW_side[qp] * vel_penalty *
                              vel_phi_side[i][qp] * vel_phi_side[j][qp];
                  if (dim == 3)
                    Kww(i,j) += JxW_side[qp] * vel_penalty *
                                vel_phi_side[i][qp] * vel_phi_side[j][qp];
                } // end of the inner dof (j) loop

            } // end - if (request_jacobian && c.elem_solution_derivative)

        } // end of the outer dof (i) loop
    } // end of the quadrature point (qp) loop

  // Pin p = p_value at p_point
  libMesh::Point p_point(0.,0.);
  if (c.elem->contains_point(p_point))
    {
      // The pressure penalty value.  \f$ \frac{1}{\epsilon} \f$
      const libMesh::Real p_penalty = 1.e9;

      libMesh::DenseSubMatrix<Number> &Kpp = *c.elem_subjacobians[p_var][p_var]; // Rp,p
      libMesh::DenseSubVector<Number> &Fp = *c.elem_subresiduals[p_var]; // Rp

      // The number of local degrees of freedom in p variable.
      const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();

      libMesh::Number p = c.point_value(p_var, p_point);
      libMesh::Number p_value = 0.;

      libMesh::FEType fe_type = c.element_fe_var[p_var]->get_fe_type();
      libMesh::Point point_loc_in_masterelem = libMesh::FEInterface::inverse_map(dim, fe_type, c.elem, zero);

      std::vector<libMesh::Real> point_p_phi(n_p_dofs);
      for (unsigned int i=0; i != n_p_dofs; i++)
          point_p_phi[i] = libMesh::FEInterface::shape(dim, fe_type, c.elem, i, point_loc_in_masterelem);

      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) += p_penalty * (p - p_value) * point_p_phi[i];
          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert (c.elem_solution_derivative == 1.0);

              for (unsigned int j=0; j != n_p_dofs; j++)
                Kpp(i,j) += p_penalty * point_p_phi[i] * point_p_phi[j];
            } // end - if (request_jacobian && c.elem_solution_derivative)
        } // end of the outer dof (i) loop
    } // end - if p_point is inside element

  return request_jacobian;
}

bool GRINS::LowMachNumberNavierStokesSystem::mass_residual(
						bool request_jacobian,
						libMesh::DiffContext& context )
{
  // TODO: account for rho factor in mass matrix
  return request_jacobian;
}

libMesh::Point GRINS::LowMachNumberNavierStokesSystem::forcing(const libMesh::Point &pt_xyz)
{
  // TODO: add forcing options
  return libMesh::Point(0.,0.,0.);
}
