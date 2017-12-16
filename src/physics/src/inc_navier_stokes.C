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
#include "grins/inc_navier_stokes.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/generic_ic_handler.h"
#include "grins/postprocessed_quantities.h"
#include "grins/inc_nav_stokes_macro.h"

// libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{

  template<class Mu>
  IncompressibleNavierStokes<Mu>::IncompressibleNavierStokes(const std::string& physics_name, const GetPot& input )
    : IncompressibleNavierStokesBase<Mu>(physics_name,
                                         PhysicsNaming::incompressible_navier_stokes(), /* "core" Physics name */
                                         input),
    _p_pinning(input,physics_name),
    _mu_index(0)
  {
    // This is deleted in the base class
    this->_ic_handler = new GenericICHandler( physics_name, input );

    this->_pin_pressure = input("Physics/"+PhysicsNaming::incompressible_navier_stokes()+"/pin_pressure", false );

    if( this->_flow_vars.dim() < 2 )
      libmesh_error_msg("ERROR: IncompressibleNavierStokes only valid for two or three dimensions! Make sure you have at least two components in your Velocity type variable.");
  }

  template<class Mu>
  void IncompressibleNavierStokes<Mu>::auxiliary_init( MultiphysicsSystem& system )
  {
    if( _pin_pressure )
      _p_pinning.check_pin_location(system.get_mesh());
  }

  template<class Mu>
  void IncompressibleNavierStokes<Mu>::register_postprocessing_vars( const GetPot& input,
                                                                     PostProcessedQuantities<libMesh::Real>& postprocessing )
  {
    std::string section = "Physics/"+PhysicsNaming::incompressible_navier_stokes()+"/output_vars";

    if( input.have_variable(section) )
      {
        unsigned int n_vars = input.vector_variable_size(section);

        for( unsigned int v = 0; v < n_vars; v++ )
          {
            std::string name = input(section,"DIE!",v);

            if( name == std::string("mu") )
              {
                this->_mu_index = postprocessing.register_quantity( name );
              }
            else
              {
                std::cerr << "Error: Invalid output_vars value for "+PhysicsNaming::incompressible_navier_stokes() << std::endl
                          << "       Found " << name << std::endl
                          << "       Acceptable values are: mu" << std::endl;
                libmesh_error();
              }
          }
      }

    return;
  }

  template<class Mu>
  void IncompressibleNavierStokes<Mu>::element_time_derivative
  ( bool compute_jacobian,
    AssemblyContext & context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();
    const unsigned int n_p_dofs = context.get_dof_indices(this->_press_var.p()).size();

    // Check number of dofs is same for this->_flow_vars.u(), v_var and w_var.
    libmesh_assert (n_u_dofs == context.get_dof_indices(this->_flow_vars.v()).size());

    if (this->_flow_vars.dim() == 3)
      libmesh_assert (n_u_dofs == context.get_dof_indices(this->_flow_vars.w()).size());

    // We get some references to cell-specific data that
    // will be used to assemble the linear system.

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // The velocity shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    // The velocity shape function gradients (in global coords.)
    // at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_flow_vars.u())->get_dphi();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& p_phi =
      context.get_element_fe(this->_press_var.p())->get_phi();

    const std::vector<libMesh::Point>& u_qpoint =
      context.get_element_fe(this->_flow_vars.u())->get_xyz();

    // The subvectors and submatrices we need to fill:
    //
    // K_{\alpha \beta} = R_{\alpha},{\beta} = \partial{ R_{\alpha} } / \partial{ {\beta} } (where R denotes residual)
    // e.g., for \alpha = v and \beta = u we get: K{vu} = R_{v},{u}
    // Note that Kpu, Kpv, Kpw and Fp comes as constraint.

    libMesh::DenseSubMatrix<libMesh::Number> &Kuu = context.get_elem_jacobian(this->_flow_vars.u(), this->_flow_vars.u()); // R_{u},{u}
    libMesh::DenseSubMatrix<libMesh::Number> &Kuv = context.get_elem_jacobian(this->_flow_vars.u(), this->_flow_vars.v()); // R_{u},{v}
    libMesh::DenseSubMatrix<libMesh::Number>* Kuw = NULL;

    libMesh::DenseSubMatrix<libMesh::Number> &Kvu = context.get_elem_jacobian(this->_flow_vars.v(), this->_flow_vars.u()); // R_{v},{u}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvv = context.get_elem_jacobian(this->_flow_vars.v(), this->_flow_vars.v()); // R_{v},{v}
    libMesh::DenseSubMatrix<libMesh::Number>* Kvw = NULL;

    libMesh::DenseSubMatrix<libMesh::Number>* Kwu = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kwv = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kww = NULL;

    libMesh::DenseSubMatrix<libMesh::Number> &Kup = context.get_elem_jacobian(this->_flow_vars.u(), this->_press_var.p()); // R_{u},{p}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvp = context.get_elem_jacobian(this->_flow_vars.v(), this->_press_var.p()); // R_{v},{p}
    libMesh::DenseSubMatrix<libMesh::Number>* Kwp = NULL;

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u()); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v()); // R_{v}
    libMesh::DenseSubVector<libMesh::Number>* Fw = NULL;

    if( this->_flow_vars.dim() == 3 )
      {
        Kuw = &context.get_elem_jacobian(this->_flow_vars.u(), this->_flow_vars.w()); // R_{u},{w}
        Kvw = &context.get_elem_jacobian(this->_flow_vars.v(), this->_flow_vars.w()); // R_{v},{w}
        Kwu = &context.get_elem_jacobian(this->_flow_vars.w(), this->_flow_vars.u()); // R_{w},{u};
        Kwv = &context.get_elem_jacobian(this->_flow_vars.w(), this->_flow_vars.v()); // R_{w},{v};
        Kww = &context.get_elem_jacobian(this->_flow_vars.w(), this->_flow_vars.w()); // R_{w},{w}
        Kwp = &context.get_elem_jacobian(this->_flow_vars.w(), this->_press_var.p()); // R_{w},{p}
        Fw  = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}
      }

    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // Compute the solution & its gradient at the old Newton iterate.
        libMesh::Number p, u, v;
        p = context.interior_value(this->_press_var.p(), qp);
        u = context.interior_value(this->_flow_vars.u(), qp);
        v = context.interior_value(this->_flow_vars.v(), qp);

        libMesh::Gradient grad_u, grad_v, grad_w;
        grad_u = context.interior_gradient(this->_flow_vars.u(), qp);
        grad_v = context.interior_gradient(this->_flow_vars.v(), qp);
        if (this->_flow_vars.dim() == 3)
          grad_w = context.interior_gradient(this->_flow_vars.w(), qp);

        libMesh::NumberVectorValue U(u,v);
        if (this->_flow_vars.dim() == 3)
          U(2) = context.interior_value(this->_flow_vars.w(), qp); // w

        const libMesh::Number  grad_u_x = grad_u(0);
        const libMesh::Number  grad_u_y = grad_u(1);
        const libMesh::Number  grad_u_z = (this->_flow_vars.dim() == 3)?grad_u(2):0;
        const libMesh::Number  grad_v_x = grad_v(0);
        const libMesh::Number  grad_v_y = grad_v(1);
        const libMesh::Number  grad_v_z = (this->_flow_vars.dim() == 3)?grad_v(2):0;
        const libMesh::Number  grad_w_x = (this->_flow_vars.dim() == 3)?grad_w(0):0;
        const libMesh::Number  grad_w_y = (this->_flow_vars.dim() == 3)?grad_w(1):0;
        const libMesh::Number  grad_w_z = (this->_flow_vars.dim() == 3)?grad_w(2):0;

        const libMesh::Number r = u_qpoint[qp](0);

        libMesh::Real jac = JxW[qp];

        // Compute the viscosity at this qp
        libMesh::Real _mu_qp = this->_mu(context, qp);

        if(Physics::is_axisymmetric())
          {
            jac *= r;
          }

        // First, an i-loop over the velocity degrees of freedom.
        // We know that n_u_dofs == n_v_dofs so we can compute contributions
        // for both at the same time.
        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) += jac *
              (-this->_rho*u_phi[i][qp]*(U*grad_u)        // convection term
               +p*u_gradphi[i][qp](0)              // pressure term
               -_mu_qp*(u_gradphi[i][qp]*grad_u) ); // diffusion term

            /*! \todo Would it be better to put this in its own DoF loop and do the if check once?*/
            if(Physics::is_axisymmetric())
              {
                Fu(i) += u_phi[i][qp]*( p/r - _mu_qp*U(0)/(r*r) )*jac;
              }

            Fv(i) += jac *
              (-this->_rho*u_phi[i][qp]*(U*grad_v)        // convection term
               +p*u_gradphi[i][qp](1)              // pressure term
               -_mu_qp*(u_gradphi[i][qp]*grad_v) ); // diffusion term

            if (this->_flow_vars.dim() == 3)
              {
                (*Fw)(i) += jac *
                  (-this->_rho*u_phi[i][qp]*(U*grad_w)        // convection term
                   +p*u_gradphi[i][qp](2)              // pressure term
                   -_mu_qp*(u_gradphi[i][qp]*grad_w) ); // diffusion term
              }

            if (compute_jacobian)
              {
                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    // TODO: precompute some terms like:
                    //   (Uvec*u_gradphi[j][qp]),
                    //   u_phi[i][qp]*u_phi[j][qp],
                    //   (u_gradphi[i][qp]*u_gradphi[j][qp])

                    Kuu(i,j) += jac * context.get_elem_solution_derivative() *
                      (-this->_rho*u_phi[i][qp]*(U*u_gradphi[j][qp])       // convection term
                       -this->_rho*u_phi[i][qp]*grad_u_x*u_phi[j][qp]             // convection term
                       -_mu_qp*(u_gradphi[i][qp]*u_gradphi[j][qp])); // diffusion term


                    if(Physics::is_axisymmetric())
                      {
                        Kuu(i,j) -= u_phi[i][qp]*_mu_qp*u_phi[j][qp]/(r*r)*jac * context.get_elem_solution_derivative();
                      }

                    Kuv(i,j) += jac * context.get_elem_solution_derivative() *
                      (-this->_rho*u_phi[i][qp]*grad_u_y*u_phi[j][qp]);           // convection term

                    Kvv(i,j) += jac * context.get_elem_solution_derivative() *
                      (-this->_rho*u_phi[i][qp]*(U*u_gradphi[j][qp])       // convection term
                       -this->_rho*u_phi[i][qp]*grad_v_y*u_phi[j][qp]             // convection term
                       -_mu_qp*(u_gradphi[i][qp]*u_gradphi[j][qp])); // diffusion term

                    Kvu(i,j) += jac * context.get_elem_solution_derivative() *
                      (-this->_rho*u_phi[i][qp]*grad_v_x*u_phi[j][qp]);           // convection term

                    if (this->_flow_vars.dim() == 3)
                      {
                        (*Kuw)(i,j) += jac * context.get_elem_solution_derivative() *
                          (-this->_rho*u_phi[i][qp]*grad_u_z*u_phi[j][qp]);           // convection term

                        (*Kvw)(i,j) += jac * context.get_elem_solution_derivative() *
                          (-this->_rho*u_phi[i][qp]*grad_v_z*u_phi[j][qp]);           // convection term

                        (*Kww)(i,j) += jac * context.get_elem_solution_derivative() *
                          (-this->_rho*u_phi[i][qp]*(U*u_gradphi[j][qp])       // convection term
                           -this->_rho*u_phi[i][qp]*grad_w_z*u_phi[j][qp]             // convection term
                           -_mu_qp*(u_gradphi[i][qp]*u_gradphi[j][qp])); // diffusion term
                        (*Kwu)(i,j) += jac * context.get_elem_solution_derivative() *
                          (-this->_rho*u_phi[i][qp]*grad_w_x*u_phi[j][qp]);           // convection term
                        (*Kwv)(i,j) += jac * context.get_elem_solution_derivative() *
                          (-this->_rho*u_phi[i][qp]*grad_w_y*u_phi[j][qp]);           // convection term
                      }
                  } // end of the inner dof (j) loop

                // Matrix contributions for the up, vp and wp couplings
                for (unsigned int j=0; j != n_p_dofs; j++)
                  {
                    Kup(i,j) += u_gradphi[i][qp](0)*p_phi[j][qp]*jac * context.get_elem_solution_derivative();
                    Kvp(i,j) += u_gradphi[i][qp](1)*p_phi[j][qp]*jac * context.get_elem_solution_derivative();

                    if (this->_flow_vars.dim() == 3)
                      {
                        (*Kwp)(i,j) += u_gradphi[i][qp](2)*p_phi[j][qp]*jac * context.get_elem_solution_derivative();
                      }

                    if(Physics::is_axisymmetric())
                      {
                        Kup(i,j) += u_phi[i][qp]*p_phi[j][qp]/r*jac * context.get_elem_solution_derivative();
                      }

                  } // end of the inner dof (j) loop



              } // end - if (compute_jacobian)

          } // end of the outer dof (i) loop
      } // end of the quadrature point (qp) loop
  }

  template<class Mu>
  void IncompressibleNavierStokes<Mu>::element_constraint
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();
    const unsigned int n_p_dofs = context.get_dof_indices(this->_press_var.p()).size();

    // We get some references to cell-specific data that
    // will be used to assemble the linear system.

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // The velocity shape function gradients (in global coords.)
    // at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_flow_vars.u())->get_dphi();

    // The velocity shape function gradients (in global coords.)
    // at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& p_phi =
      context.get_element_fe(this->_press_var.p())->get_phi();

    const std::vector<libMesh::Point>& u_qpoint =
      context.get_element_fe(this->_flow_vars.u())->get_xyz();

    // The subvectors and submatrices we need to fill:
    //
    // Kpu, Kpv, Kpw, Fp

    libMesh::DenseSubMatrix<libMesh::Number> &Kpu = context.get_elem_jacobian(this->_press_var.p(), this->_flow_vars.u()); // R_{p},{u}
    libMesh::DenseSubMatrix<libMesh::Number> &Kpv = context.get_elem_jacobian(this->_press_var.p(), this->_flow_vars.v()); // R_{p},{v}
    libMesh::DenseSubMatrix<libMesh::Number>* Kpw = NULL;

    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(this->_press_var.p()); // R_{p}


    if( this->_flow_vars.dim() == 3 )
      {
        Kpw = &context.get_elem_jacobian(this->_press_var.p(), this->_flow_vars.w()); // R_{p},{w}
      }

    // Add the constraint given by the continuity equation.
    unsigned int n_qpoints = context.get_element_qrule().n_points();
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // Compute the velocity gradient at the old Newton iterate.
        libMesh::Gradient grad_u, grad_v, grad_w;
        grad_u = context.interior_gradient(this->_flow_vars.u(), qp);
        grad_v = context.interior_gradient(this->_flow_vars.v(), qp);
        if (this->_flow_vars.dim() == 3)
          grad_w = context.interior_gradient(this->_flow_vars.w(), qp);

        libMesh::Number divU = grad_u(0) + grad_v(1);
        if (this->_flow_vars.dim() == 3)
          divU += grad_w(2);

        const libMesh::Number r = u_qpoint[qp](0);

        libMesh::Real jac = JxW[qp];

        if(Physics::is_axisymmetric())
          {
            libMesh::Number u = context.interior_value( this->_flow_vars.u(), qp );
            divU += u/r;
            jac *= r;
          }

        // Now a loop over the pressure degrees of freedom.  This
        // computes the contributions of the continuity equation.
        for (unsigned int i=0; i != n_p_dofs; i++)
          {
            Fp(i) += p_phi[i][qp]*divU*jac;

            if (compute_jacobian)
              {
                libmesh_assert_equal_to (context.get_elem_solution_derivative(), 1.0);

                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    Kpu(i,j) += p_phi[i][qp]*u_gradphi[j][qp](0)*jac * context.get_elem_solution_derivative();
                    Kpv(i,j) += p_phi[i][qp]*u_gradphi[j][qp](1)*jac * context.get_elem_solution_derivative();
                    if (this->_flow_vars.dim() == 3)
                      (*Kpw)(i,j) += p_phi[i][qp]*u_gradphi[j][qp](2)*jac * context.get_elem_solution_derivative();

                    if(Physics::is_axisymmetric())
                      {
                        Kpu(i,j) += p_phi[i][qp]*u_phi[j][qp]/r*jac * context.get_elem_solution_derivative();
                      }
                  } // end of the inner dof (j) loop

              } // end - if (compute_jacobian)

          } // end of the outer dof (i) loop
      } // end of the quadrature point (qp) loop

    // Pin p = p_value at p_point
    if( _pin_pressure )
      {
        _p_pinning.pin_value( context, compute_jacobian, this->_press_var.p() );
      }
  }

  template<class Mu>
  void IncompressibleNavierStokes<Mu>::mass_residual
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // Element Jacobian * quadrature weights for interior integration
    // We assume the same for each flow variable
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // The shape functions at interior quadrature points.
    // We assume the same for each flow variable
    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    const std::vector<libMesh::Point>& u_qpoint =
      context.get_element_fe(this->_flow_vars.u())->get_xyz();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<libMesh::Real> &F_u = context.get_elem_residual(this->_flow_vars.u());
    libMesh::DenseSubVector<libMesh::Real> &F_v = context.get_elem_residual(this->_flow_vars.v());
    libMesh::DenseSubVector<libMesh::Real>* F_w = NULL;

    libMesh::DenseSubMatrix<libMesh::Real> &M_uu = context.get_elem_jacobian(this->_flow_vars.u(), this->_flow_vars.u());
    libMesh::DenseSubMatrix<libMesh::Real> &M_vv = context.get_elem_jacobian(this->_flow_vars.v(), this->_flow_vars.v());
    libMesh::DenseSubMatrix<libMesh::Real>* M_ww = NULL;


    if( this->_flow_vars.dim() == 3 )
      {
        F_w  = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}
        M_ww = &context.get_elem_jacobian(this->_flow_vars.w(), this->_flow_vars.w());
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
        // For the mass residual, we need to be a little careful.
        // The time integrator is handling the time-discretization
        // for us so we need to supply M(u_fixed)*u' for the residual.
        // u_fixed will be given by the fixed_interior_value function
        // while u' will be given by the interior_rate function.
        libMesh::Real u_dot, v_dot, w_dot = 0.0;
        context.interior_rate(this->_flow_vars.u(), qp, u_dot);
        context.interior_rate(this->_flow_vars.v(), qp, v_dot);

        if(this->_flow_vars.dim() == 3 )
          context.interior_rate(this->_flow_vars.w(), qp, w_dot);

        const libMesh::Number r = u_qpoint[qp](0);

        libMesh::Real jac = JxW[qp];

        if(Physics::is_axisymmetric())
          {
            jac *= r;
          }

        for (unsigned int i = 0; i != n_u_dofs; ++i)
          {
            F_u(i) -= this->_rho*u_dot*u_phi[i][qp]*jac;
            F_v(i) -= this->_rho*v_dot*u_phi[i][qp]*jac;

            if( this->_flow_vars.dim() == 3 )
              (*F_w)(i) -= this->_rho*w_dot*u_phi[i][qp]*jac;

            if( compute_jacobian )
              {
                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    // Assuming rho is constant w.r.t. u, v, w
                    // and T (if Boussinesq added).
                    libMesh::Real value = this->_rho*u_phi[i][qp]*u_phi[j][qp]*jac * context.get_elem_solution_rate_derivative();

                    M_uu(i,j) -= value;
                    M_vv(i,j) -= value;

                    if( this->_flow_vars.dim() == 3)
                      {
                        (*M_ww)(i,j) -= value;
                      }

                  } // End dof loop
              } // End Jacobian check
          } // End dof loop
      } // End quadrature loop

    return;
  }

  template<class Mu>
  void IncompressibleNavierStokes<Mu>::compute_postprocessed_quantity( unsigned int quantity_index,
                                                                       const AssemblyContext& context,
                                                                       const libMesh::Point& point,
                                                                       libMesh::Real& value )
  {
    if( quantity_index == this->_mu_index )
      {
        value = this->_mu(point, context.get_time());
      }

    return;
  }

} // namespace GRINS

// Instantiate
INSTANTIATE_INC_NS_SUBCLASS(IncompressibleNavierStokes);
