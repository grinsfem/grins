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
#include "grins/low_mach_navier_stokes.h"

// GRINS
#include "grins_config.h"
#include "grins/assembly_context.h"
#include "grins/constant_viscosity.h"
#include "grins/constant_specific_heat.h"
#include "grins/constant_conductivity.h"
#include "grins/generic_ic_handler.h"
#include "grins/postprocessed_quantities.h"

// libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{

  template<class Mu, class SH, class TC>
  LowMachNavierStokes<Mu,SH,TC>::LowMachNavierStokes(const std::string& physics_name, const GetPot& input)
    : LowMachNavierStokesBase<Mu,SH,TC>(physics_name,PhysicsNaming::low_mach_navier_stokes(),input),
    _p_pinning(input,physics_name),
    _rho_index(0) // Initialize to zero
  {
    this->_ic_handler = new GenericICHandler( physics_name, input );

    this->_pin_pressure = input("Physics/"+PhysicsNaming::low_mach_navier_stokes()+"/pin_pressure", false );

    if( this->_flow_vars.dim() < 2 )
      libmesh_error_msg("ERROR: LowMachNavierStokes only valid for two or three dimensions! Make sure you have at least two components in your Velocity type variable.");
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokes<Mu,SH,TC>::auxiliary_init( MultiphysicsSystem& system )
  {
    if( _pin_pressure )
      _p_pinning.check_pin_location(system.get_mesh());
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokes<Mu,SH,TC>::register_postprocessing_vars( const GetPot& input,
                                                                    PostProcessedQuantities<libMesh::Real>& postprocessing )
  {
    std::string section = "Physics/"+PhysicsNaming::low_mach_navier_stokes()+"/output_vars";

    if( input.have_variable(section) )
      {
        unsigned int n_vars = input.vector_variable_size(section);

        for( unsigned int v = 0; v < n_vars; v++ )
          {
            std::string name = input(section,"DIE!",v);

            if( name == std::string("rho") )
              {
                this->_rho_index = postprocessing.register_quantity( name );
              }
            else
              {
                std::cerr << "Error: Invalue output_vars value for "+PhysicsNaming::low_mach_navier_stokes() << std::endl
                          << "       Found " << name << std::endl
                          << "       Acceptable values are: rho" << std::endl;
                libmesh_error();
              }
          }
      }
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokes<Mu,SH,TC>::init_context( AssemblyContext & context )
  {
    // First call base class
    LowMachNavierStokesBase<Mu,SH,TC>::init_context(context);

    // We also need the side shape functions, etc.
    context.get_side_fe(this->_flow_vars.u())->get_JxW();
    context.get_side_fe(this->_flow_vars.u())->get_phi();
    context.get_side_fe(this->_flow_vars.u())->get_dphi();
    context.get_side_fe(this->_flow_vars.u())->get_xyz();

    context.get_side_fe(this->_temp_vars.T())->get_JxW();
    context.get_side_fe(this->_temp_vars.T())->get_phi();
    context.get_side_fe(this->_temp_vars.T())->get_dphi();
    context.get_side_fe(this->_temp_vars.T())->get_xyz();
  }


  template<class Mu, class SH, class TC>
  void LowMachNavierStokes<Mu,SH,TC>::element_time_derivative
  ( bool compute_jacobian, AssemblyContext & context )
  {
    const CachedValues & cache = context.get_cached_values();

    this->assemble_mass_time_deriv( compute_jacobian, context, cache );
    this->assemble_momentum_time_deriv( compute_jacobian, context, cache );
    this->assemble_energy_time_deriv( compute_jacobian, context, cache );

    if( this->_enable_thermo_press_calc )
      this->assemble_thermo_press_elem_time_deriv( compute_jacobian, context );
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokes<Mu,SH,TC>::element_constraint
  ( bool compute_jacobian,
    AssemblyContext & context )
  {
    // Pin p = p_value at p_point
    if( this->_pin_pressure )
      this->_p_pinning.pin_value( context, compute_jacobian, this->_press_var.p());
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokes<Mu,SH,TC>::mass_residual
  ( bool compute_jacobian, AssemblyContext & context )
  {
    this->assemble_continuity_mass_residual( compute_jacobian, context );

    this->assemble_momentum_mass_residual( compute_jacobian, context );

    this->assemble_energy_mass_residual( compute_jacobian, context );

    if( this->_enable_thermo_press_calc )
      this->assemble_thermo_press_mass_residual( compute_jacobian, context );
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokes<Mu,SH,TC>::assemble_mass_time_deriv( bool compute_jacobian,
                                                                AssemblyContext & context,
                                                                const CachedValues & cache )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_p_dofs = context.get_dof_indices(this->_press_var.p()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& p_phi =
      context.get_element_fe(this->_press_var.p())->get_phi();



    // The velocity shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    // The velocity shape function gradients (in global coords.)
    // at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_flow_vars.u())->get_dphi();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_temp_vars.T())->get_phi();

    // The temperature shape function gradients (in global coords.)
    // at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(this->_temp_vars.T())->get_dphi();



    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(this->_press_var.p()); // R_{p}

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    // The number of local degrees of freedom in each variable.
    const unsigned int n_t_dofs = context.get_dof_indices(this->_temp_vars.T()).size();

    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();

    // Check number of dofs is same for _flow_vars.u(), v_var and w_var.
    libmesh_assert (n_u_dofs == context.get_dof_indices(this->_flow_vars.v()).size());

    if (this->_flow_vars.dim() == 3)
      libmesh_assert (n_u_dofs == context.get_dof_indices(this->_flow_vars.w()).size());

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Number u, v, T;
        u = cache.get_cached_values(Cache::X_VELOCITY)[qp];
        v = cache.get_cached_values(Cache::Y_VELOCITY)[qp];

        T = cache.get_cached_values(Cache::TEMPERATURE)[qp];

        libMesh::Gradient grad_u = cache.get_cached_gradient_values(Cache::X_VELOCITY_GRAD)[qp];
        libMesh::Gradient grad_v = cache.get_cached_gradient_values(Cache::Y_VELOCITY_GRAD)[qp];

        libMesh::Gradient grad_T = cache.get_cached_gradient_values(Cache::TEMPERATURE_GRAD)[qp];

        libMesh::NumberVectorValue U(u,v);
        if (this->_flow_vars.dim() == 3)
          U(2) = cache.get_cached_values(Cache::Z_VELOCITY)[qp]; // w

        libMesh::Number divU = grad_u(0) + grad_v(1);
        if (this->_flow_vars.dim() == 3)
          {
            libMesh::Gradient grad_w = cache.get_cached_gradient_values(Cache::Z_VELOCITY_GRAD)[qp];
            divU += grad_w(2);
          }


        libMesh::DenseSubMatrix<libMesh::Number> &KPu = context.get_elem_jacobian(this->_press_var.p(), this->_flow_vars.u());
        libMesh::DenseSubMatrix<libMesh::Number> &KPv = context.get_elem_jacobian(this->_press_var.p(), this->_flow_vars.v());
        libMesh::DenseSubMatrix<libMesh::Number> &KPT = context.get_elem_jacobian(this->_press_var.p(), this->_temp_vars.T());

        libMesh::DenseSubMatrix<libMesh::Number>* KPw = NULL;

        if( this->_flow_vars.dim() == 3 )
          {
            KPw = &context.get_elem_jacobian(this->_press_var.p(), this->_flow_vars.w());
          }

        // Now a loop over the pressure degrees of freedom.  This
        // computes the contributions of the continuity equation.
        for (unsigned int i=0; i != n_p_dofs; i++)
          {
            Fp(i) += (-U*grad_T/T + divU)*p_phi[i][qp]*JxW[qp];


            if (compute_jacobian)
              {

                for (unsigned int j=0; j!=n_u_dofs; j++)
                  {
                    KPu(i,j) += JxW[qp]*(
                                         +u_gradphi[j][qp](0)*p_phi[i][qp]
                                         -u_phi[j][qp]*p_phi[i][qp]*grad_T(0)/T
                                         );

                    KPv(i,j) += JxW[qp]*(
                                         +u_gradphi[j][qp](1)*p_phi[i][qp]
                                         -u_phi[j][qp]*p_phi[i][qp]*grad_T(1)/T
                                         );

                    if (this->_flow_vars.dim() == 3)
                      {
                        (*KPw)(i,j) += JxW[qp]*(
                                                +u_gradphi[j][qp](2)*p_phi[i][qp]
                                                -u_phi[j][qp]*p_phi[i][qp]*grad_T(2)/T
                                                );
                      }

                  }

                for (unsigned int j=0; j!=n_t_dofs; j++)
                  {
                    KPT(i,j) += JxW[qp]*(
                                         -T_gradphi[j][qp]*U*p_phi[i][qp]/T
                                         +U*p_phi[i][qp]*grad_T*T_phi[j][qp])/(T*T);
                  }
              } // end if compute_jacobian
          } // end p_dofs loop
      } // end qp loop

    return;
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokes<Mu,SH,TC>::assemble_momentum_time_deriv( bool compute_jacobian,
                                                                    AssemblyContext & context,
                                                                    const CachedValues & cache )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();
    const unsigned int n_p_dofs = context.get_dof_indices(this->_press_var.p()).size();
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();

    // Check number of dofs is same for _flow_vars.u(), v_var and w_var.
    libmesh_assert (n_u_dofs == context.get_dof_indices(this->_flow_vars.v()).size());

    if (this->_flow_vars.dim() == 3)
      libmesh_assert (n_u_dofs == context.get_dof_indices(this->_flow_vars.w()).size());

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();
    const std::vector<std::vector<libMesh::Real> >& p_phi =
      context.get_element_fe(this->_press_var.p())->get_phi();
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_temp_vars.T())->get_phi();

    // The velocity shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_flow_vars.u())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u()); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v()); // R_{v}
    libMesh::DenseSubVector<libMesh::Real>* Fw = NULL;

    if( this->_flow_vars.dim() == 3 )
      {
        Fw  = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}
      }

    unsigned int n_qpoints = context.get_element_qrule().n_points();
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Number u, v, p, p0, T;
        u = cache.get_cached_values(Cache::X_VELOCITY)[qp];
        v = cache.get_cached_values(Cache::Y_VELOCITY)[qp];

        T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
        p = cache.get_cached_values(Cache::PRESSURE)[qp];
        p0 = cache.get_cached_values(Cache::THERMO_PRESSURE)[qp];

        libMesh::Gradient grad_u = cache.get_cached_gradient_values(Cache::X_VELOCITY_GRAD)[qp];
        libMesh::Gradient grad_v = cache.get_cached_gradient_values(Cache::Y_VELOCITY_GRAD)[qp];

        libMesh::Gradient grad_w;
        if (this->_flow_vars.dim() == 3)
          grad_w = cache.get_cached_gradient_values(Cache::Z_VELOCITY_GRAD)[qp];

        libMesh::NumberVectorValue grad_uT( grad_u(0), grad_v(0) );
        libMesh::NumberVectorValue grad_vT( grad_u(1), grad_v(1) );
        libMesh::NumberVectorValue grad_wT;
        if( this->_flow_vars.dim() == 3 )
          {
            grad_uT(2) = grad_w(0);
            grad_vT(2) = grad_w(1);
            grad_wT = libMesh::NumberVectorValue( grad_u(2), grad_v(2), grad_w(2) );
          }

        libMesh::NumberVectorValue U(u,v);
        if (this->_flow_vars.dim() == 3)
          U(2) = cache.get_cached_values(Cache::Z_VELOCITY)[qp]; // w

        libMesh::Number divU = grad_u(0) + grad_v(1);
        if (this->_flow_vars.dim() == 3)
          divU += grad_w(2);

        libMesh::Number rho = this->rho( T, p0 );
        libMesh::Number d_rho = this->d_rho_dT( T, p0 );
        libMesh::Number d_mu = this->_mu.deriv(T);

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

        libMesh::DenseSubMatrix<libMesh::Number> &KuT = context.get_elem_jacobian(this->_flow_vars.u(), this->_temp_vars.T()); // R_{u},{p}
        libMesh::DenseSubMatrix<libMesh::Number> &KvT = context.get_elem_jacobian(this->_flow_vars.v(), this->_temp_vars.T()); // R_{v},{p}
        libMesh::DenseSubMatrix<libMesh::Number>* KwT = NULL;

        if( this->_flow_vars.dim() == 3 )
          {
            Kuw = &context.get_elem_jacobian(this->_flow_vars.u(), this->_flow_vars.w()); // R_{u},{w}
            Kvw = &context.get_elem_jacobian(this->_flow_vars.v(), this->_flow_vars.w()); // R_{v},{w}
            Kwu = &context.get_elem_jacobian(this->_flow_vars.w(), this->_flow_vars.u()); // R_{w},{u};
            Kwv = &context.get_elem_jacobian(this->_flow_vars.w(), this->_flow_vars.v()); // R_{w},{v};
            Kww = &context.get_elem_jacobian(this->_flow_vars.w(), this->_flow_vars.w()); // R_{w},{w}
            Kwp = &context.get_elem_jacobian(this->_flow_vars.w(), this->_press_var.p()); // R_{w},{p}
            KwT = &context.get_elem_jacobian(this->_flow_vars.w(), this->_temp_vars.T()); // R_{w},{T}
          }

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
            if (this->_flow_vars.dim() == 3)
              {
                (*Fw)(i) += ( -rho*U*grad_w*u_phi[i][qp]                 // convection term
                              + p*u_gradphi[i][qp](2)                           // pressure term
                              - this->_mu(T)*(u_gradphi[i][qp]*grad_w + u_gradphi[i][qp]*grad_wT
                                              - 2.0/3.0*divU*u_gradphi[i][qp](2) )    // diffusion term
                              + rho*this->_g(2)*u_phi[i][qp]                 // hydrostatic term
                              )*JxW[qp];
              }

            if (compute_jacobian && context.get_elem_solution_derivative())
              {
                libmesh_assert (context.get_elem_solution_derivative() == 1.0);

                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    //precompute repeated terms
                    libMesh::Number r0 = rho*U*u_phi[i][qp]*u_gradphi[j][qp];
                    libMesh::Number r1 = u_gradphi[i][qp]*u_gradphi[j][qp];
                    libMesh::Number r2 = rho*u_phi[i][qp]*u_phi[j][qp];


                    Kuu(i,j) += JxW[qp]*(
                                         -r0
                                         //-rho*U*u_gradphi[j][qp]*u_phi[i][qp]
                                         -r2*grad_u(0)
                                         //-rho*u_phi[i][qp]*grad_u(0)*u_phi[j][qp]
                                         -this->_mu(T)*(
                                                        r1
                                                        //u_gradphi[i][qp]*u_gradphi[j][qp]
                                                        + u_gradphi[i][qp](0)*u_gradphi[j][qp](0) // transpose
                                                        - 2.0/3.0*u_gradphi[i][qp](0)*u_gradphi[j][qp](0)
                                                        ));
                    Kvv(i,j) += JxW[qp]*(
                                         -r0
                                         //-rho*U*u_gradphi[j][qp]*u_phi[i][qp]
                                         -r2*grad_v(1)
                                         //-rho*u_phi[i][qp]*grad_v(1)*u_phi[j][qp]
                                         -this->_mu(T)*(
                                                        r1
                                                        //u_gradphi[i][qp]*u_gradphi[j][qp]
                                                        + u_gradphi[i][qp](1)*u_gradphi[j][qp](1) // transpose
                                                        - 2.0/3.0*u_gradphi[i][qp](1)*u_gradphi[j][qp](1)
                                                        ));

                    Kuv(i,j) += JxW[qp]*(
                                         +2.0/3.0*this->_mu(T)*u_gradphi[i][qp](0)*u_gradphi[j][qp](1)
                                         -this->_mu(T)*u_gradphi[i][qp](1)*u_gradphi[j][qp](0)
                                         -r2*grad_u(1)
                                         //-rho*u_phi[i][qp]*u_phi[j][qp]*grad_u(1));
                                         );

                    Kvu(i,j) += JxW[qp]*(
                                         +2.0/3.0*this->_mu(T)*u_gradphi[i][qp](1)*u_gradphi[j][qp](0)
                                         -this->_mu(T)*u_gradphi[i][qp](0)*u_gradphi[j][qp](1)
                                         -r2*grad_v(0)
                                         //-rho*u_phi[i][qp]*u_phi[j][qp]*grad_v(0));
                                         );



                    if (this->_flow_vars.dim() == 3)
                      {
                        (*Kuw)(i,j) += JxW[qp]*(
                                                +2.0/3.0*this->_mu(T)*u_gradphi[i][qp](0)*u_gradphi[j][qp](2)
                                                -this->_mu(T)*u_gradphi[i][qp](2)*u_gradphi[j][qp](0)
                                                -r2*grad_u(2)
                                                //-rho*u_phi[i][qp]*u_phi[j][qp]*grad_u(2)
                                                );

                        (*Kvw)(i,j) += JxW[qp]*(
                                                +2.0/3.0*this->_mu(T)*u_gradphi[i][qp](1)*u_gradphi[j][qp](2)
                                                -this->_mu(T)*u_gradphi[i][qp](2)*u_gradphi[j][qp](1)
                                                -r2*grad_v(2)
                                                //-rho*u_phi[i][qp]*u_phi[j][qp]*grad_v(2)
                                                );

                        (*Kwu)(i,j) += JxW[qp]*(
                                                +2.0/3.0*this->_mu(T)*u_gradphi[i][qp](2)*u_gradphi[j][qp](0)
                                                -this->_mu(T)*u_gradphi[i][qp](0)*u_gradphi[j][qp](2)
                                                -r2*grad_w(0)
                                                //-rho*u_phi[i][qp]*u_phi[j][qp]*grad_w(0)
                                                );

                        (*Kwv)(i,j) += JxW[qp]*(
                                                +2.0/3.0*this->_mu(T)*u_gradphi[i][qp](2)*u_gradphi[j][qp](1)
                                                -this->_mu(T)*u_gradphi[i][qp](1)*u_gradphi[j][qp](2)
                                                -r2*grad_w(1)
                                                //-rho*u_phi[i][qp]*u_phi[j][qp]*grad_w(1)
                                                );

                        (*Kww)(i,j) += JxW[qp]*(
                                                -r0
                                                //-rho*U*u_gradphi[j][qp]*u_phi[i][qp]
                                                -r2*grad_w(2)
                                                //-rho*u_phi[i][qp]*grad_w(2)*u_phi[j][qp]
                                                -this->_mu(T)*(
                                                               r1
                                                               //u_gradphi[i][qp]*u_gradphi[j][qp]
                                                               + u_gradphi[i][qp](2)*u_gradphi[j][qp](2) // transpose
                                                               - 2.0/3.0*u_gradphi[i][qp](2)*u_gradphi[j][qp](2)
                                                               ));

                      } // end if _dim==3
                  } // end of the inner dof (j) loop

                for (unsigned int j=0; j!=n_T_dofs; j++)
                  {

                    //precompute repeated term
                    libMesh:: Number r3 = d_rho*u_phi[i][qp]*T_phi[j][qp];

                    // Analytical Jacobains
                    KuT(i,j) += JxW[qp]*(
                                         -r3*U*grad_u
                                         //-d_rho*T_phi[j][qp]*U*grad_u*u_phi[i][qp]
                                         -d_mu*T_phi[j][qp]*grad_u*u_gradphi[i][qp]
                                         -d_mu*T_phi[j][qp]*grad_u(0)*u_gradphi[i][qp](0) // transpose
                                         +2.0/3.0*d_mu*T_phi[j][qp]*divU*u_gradphi[i][qp](0)
                                         +r3*this->_g(0)
                                         //+d_rho*T_phi[j][qp]*this->_g(0)*u_phi[i][qp]
                                         );

                    KvT(i,j) += JxW[qp]*(
                                         -r3*U*grad_v
                                         //-d_rho*T_phi[j][qp]*U*grad_v*u_phi[i][qp]
                                         -d_mu*T_phi[j][qp]*grad_v*u_gradphi[i][qp]
                                         -d_mu*T_phi[j][qp]*grad_v(1)*u_gradphi[i][qp](1) // transpose
                                         +2.0/3.0*d_mu*T_phi[j][qp]*divU*u_gradphi[i][qp](1)
                                         +r3*this->_g(1)
                                         //+d_rho*T_phi[j][qp]*this->_g(1)*u_phi[i][qp]
                                         );

                    if (this->_flow_vars.dim() == 3)
                      {
                        (*KwT)(i,j) += JxW[qp]*(
                                                -r3*U*grad_w
                                                //-d_rho*T_phi[j][qp]*U*grad_v*u_phi[i][qp]
                                                -d_mu*T_phi[j][qp]*grad_w*u_gradphi[i][qp]
                                                -d_mu*T_phi[j][qp]*grad_w(2)*u_gradphi[i][qp](2) // transpose
                                                +2.0/3.0*d_mu*T_phi[j][qp]*divU*u_gradphi[i][qp](2)
                                                +r3*this->_g(2)
                                                //+d_rho*T_phi[j][qp]*this->_g(2)*u_phi[i][qp]
                                                );

                      } // end if _dim==3
                  } // end T_dofs loop

                // Matrix contributions for the up, vp and wp couplings
                for (unsigned int j=0; j != n_p_dofs; j++)
                  {
                    Kup(i,j) += JxW[qp]*p_phi[j][qp]*u_gradphi[i][qp](0);
                    Kvp(i,j) += JxW[qp]*p_phi[j][qp]*u_gradphi[i][qp](1);
                    if (this->_flow_vars.dim() == 3)
                      (*Kwp)(i,j) += JxW[qp]*p_phi[j][qp]*u_gradphi[i][qp](2);
                  } // end of the inner dof (j) loop

              } // end - if (compute_jacobian && context.get_elem_solution_derivative())

          } // end of the outer dof (i) loop
      } // end of the quadrature point (qp) loop

    return;
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokes<Mu,SH,TC>::assemble_energy_time_deriv( bool compute_jacobian,
                                                                  AssemblyContext & context,
                                                                  const CachedValues & cache )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_temp_vars.T())->get_JxW();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_temp_vars.T())->get_phi();
    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    // The temperature shape functions gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(this->_temp_vars.T())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &FT = context.get_elem_residual(this->_temp_vars.T()); // R_{T}

    unsigned int n_qpoints = context.get_element_qrule().n_points();
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Number u, v, T, p0;
        u = cache.get_cached_values(Cache::X_VELOCITY)[qp];
        v = cache.get_cached_values(Cache::Y_VELOCITY)[qp];
        T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
        p0 = cache.get_cached_values(Cache::THERMO_PRESSURE)[qp];

        libMesh::Gradient grad_T = cache.get_cached_gradient_values(Cache::TEMPERATURE_GRAD)[qp];

        libMesh::NumberVectorValue U(u,v);
        if (this->_flow_vars.dim() == 3)
          U(2) = cache.get_cached_values(Cache::Z_VELOCITY)[qp]; // w

        libMesh::DenseSubMatrix<libMesh::Number> &KTu = context.get_elem_jacobian(this->_temp_vars.T(), this->_flow_vars.u()); // R_{u},{u}
        libMesh::DenseSubMatrix<libMesh::Number> &KTv = context.get_elem_jacobian(this->_temp_vars.T(), this->_flow_vars.v()); // R_{u},{u}
        libMesh::DenseSubMatrix<libMesh::Number>* KTw = NULL;

        libMesh::DenseSubMatrix<libMesh::Number> &KTT = context.get_elem_jacobian(this->_temp_vars.T(), this->_temp_vars.T()); // R_{u},{u}

        if( this->_flow_vars.dim() == 3 )
          {
            KTw = &context.get_elem_jacobian(this->_temp_vars.T(), this->_flow_vars.w()); // R_{u},{w}
          }

        libMesh::Number k = this->_k(T);
        libMesh::Number dk_dT = this->_k.deriv(T);
        libMesh::Number cp = this->_cp(T);
        libMesh::Number d_cp = this->_cp.deriv(T);
        libMesh::Number rho = this->rho( T, p0 );
        libMesh::Number d_rho = this->d_rho_dT( T, p0 );

        // Now a loop over the pressure degrees of freedom.  This
        // computes the contributions of the continuity equation.
        for (unsigned int i=0; i != n_T_dofs; i++)
          {
            FT(i) += ( -rho*cp*U*grad_T*T_phi[i][qp] // convection term
                       - k*grad_T*T_gradphi[i][qp]            // diffusion term
                       )*JxW[qp];


            if(compute_jacobian)
              {
                for (unsigned int j=0; j!=n_u_dofs; j++)
                  {
                    //pre-compute repeated term
                    libMesh::Number r0 = rho*cp*T_phi[i][qp]*u_phi[j][qp];

                    KTu(i,j) += JxW[qp]*
                      -r0*grad_T(0);
                    //-rho*cp*u_phi[j][qp]*grad_T(0)*T_phi[i][qp];

                    KTv(i,j) += JxW[qp]*
                      -r0*grad_T(1);
                    //-rho*cp*u_phi[j][qp]*grad_T(1)*T_phi[i][qp];

                    if (this->_flow_vars.dim() == 3)
                      {
                        (*KTw)(i,j) += JxW[qp]*
                          -r0*grad_T(2);
                        //-rho*cp*u_phi[j][qp]*grad_T(2)*T_phi[i][qp];
                      }

                  } // end u_dofs loop (j)


                for (unsigned int j=0; j!=n_T_dofs; j++)
                  {
                    KTT(i,j) += JxW[qp]* (
                                          -rho*(
                                                cp*U*T_phi[i][qp]*T_gradphi[j][qp]
                                                + U*grad_T*T_phi[i][qp]*d_cp*T_phi[j][qp]
                                                )
                                          -cp*U*grad_T*T_phi[i][qp]*d_rho*T_phi[j][qp]
                                          -k*T_gradphi[i][qp]*T_gradphi[j][qp]
                                          -grad_T*T_gradphi[i][qp]*dk_dT*T_phi[j][qp]
                                          );
                  } // end T_dofs loop (j)

              } // end if compute_jacobian
          } // end outer T_dofs loop (i)
      } //end qp loop

    return;
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokes<Mu,SH,TC>::assemble_continuity_mass_residual( bool /*compute_jacobian*/,
                                                                         AssemblyContext& context )
  {
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // The shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& p_phi =
      context.get_element_fe(this->_press_var.p())->get_phi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_p_dofs = context.get_dof_indices(this->_press_var.p()).size();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<libMesh::Real> &F_p = context.get_elem_residual(this->_press_var.p());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
        // For the mass residual, we need to be a little careful.
        // The time integrator is handling the time-discretization
        // for us so we need to supply M(u_fixed)*u' for the residual.
        // u_fixed will be given by the fixed_interior_value function
        // while u' will be given by the interior_rate function.
        libMesh::Real T_dot;
        context.interior_rate(this->_temp_vars.T(), qp, T_dot);

        libMesh::Real T = context.fixed_interior_value(this->_temp_vars.T(), qp);

        for (unsigned int i = 0; i != n_p_dofs; ++i)
          {
            F_p(i) -= T_dot/T*p_phi[i][qp]*JxW[qp];
          } // End DoF loop i

      } // End quadrature loop qp

    return;
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokes<Mu,SH,TC>::assemble_momentum_mass_residual( bool /*compute_jacobian*/,
                                                                       AssemblyContext& context )
  {
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // The shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u()).size();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<libMesh::Real> &F_u = context.get_elem_residual(this->_flow_vars.u());
    libMesh::DenseSubVector<libMesh::Real> &F_v = context.get_elem_residual(this->_flow_vars.v());
    libMesh::DenseSubVector<libMesh::Real>* F_w = NULL;


    if( this->_flow_vars.dim() == 3 )
      F_w  = &context.get_elem_residual(this->_flow_vars.w()); // R_{w}

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

        if( this->_flow_vars.dim() == 3 )
          context.interior_rate(this->_flow_vars.w(), qp, w_dot);

        libMesh::Real T = context.fixed_interior_value(this->_temp_vars.T(), qp);

        libMesh::Number rho = this->rho(T, this->get_p0_transient(context, qp));

        for (unsigned int i = 0; i != n_u_dofs; ++i)
          {
            F_u(i) -= rho*u_dot*u_phi[i][qp]*JxW[qp];
            F_v(i) -= rho*v_dot*u_phi[i][qp]*JxW[qp];

            if( this->_flow_vars.dim() == 3 )
              (*F_w)(i) -= rho*w_dot*u_phi[i][qp]*JxW[qp];

            /*
              if( compute_jacobian )
              {
              for (unsigned int j=0; j != n_u_dofs; j++)
              {
              // Assuming rho is constant w.r.t. u, v, w
              // and T (if Boussinesq added).
              libMesh::Real value = JxW[qp]*_rho*u_phi[i][qp]*u_phi[j][qp];

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
  void LowMachNavierStokes<Mu,SH,TC>::assemble_energy_mass_residual( bool /*compute_jacobian*/,
                                                                     AssemblyContext& context )
  {
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_temp_vars.T())->get_JxW();

    // The shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_temp_vars.T())->get_phi();

    // The number of local degrees of freedom in each variable
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<libMesh::Real> &F_T = context.get_elem_residual(this->_temp_vars.T());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
        // For the mass residual, we need to be a little careful.
        // The time integrator is handling the time-discretization
        // for us so we need to supply M(u_fixed)*u' for the residual.
        // u_fixed will be given by the fixed_interior_value function
        // while u will be given by the interior_rate function.
        libMesh::Real T_dot;
        context.interior_rate(this->_temp_vars.T(), qp, T_dot);

        libMesh::Real T = context.fixed_interior_value(this->_temp_vars.T(), qp);

        libMesh::Real cp = this->_cp(T);

        libMesh::Number rho = this->rho(T, this->get_p0_transient(context, qp));

        for (unsigned int i = 0; i != n_T_dofs; ++i)
          {
            F_T(i) -= rho*cp*T_dot*T_phi[i][qp]*JxW[qp];
          } // End DoF loop i

      } // End quadrature loop qp

    return;
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokes<Mu,SH,TC>::assemble_thermo_press_elem_time_deriv( bool /*compute_jacobian*/,
                                                                             AssemblyContext& context )
  {
    // Element Jacobian * quadrature weights for interior integration
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_temp_vars.T())->get_JxW();

    // The number of local degrees of freedom in each variable
    const unsigned int n_p0_dofs = context.get_dof_indices(this->_p0_var->p0()).size();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<libMesh::Real> &F_p0 = context.get_elem_residual(this->_p0_var->p0());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
        libMesh::Number T;
        T = context.interior_value(this->_temp_vars.T(), qp);

        libMesh::Gradient grad_u, grad_v, grad_w;
        grad_u = context.interior_gradient(this->_flow_vars.u(), qp);
        grad_v = context.interior_gradient(this->_flow_vars.v(), qp);
        if (this->_flow_vars.dim() == 3)
          grad_w = context.interior_gradient(this->_flow_vars.w(), qp);

        libMesh::Number divU = grad_u(0) + grad_v(1);
        if(this->_flow_vars.dim()==3)
          divU += grad_w(2);

        //libMesh::Number cp = this->_cp(T);
        //libMesh::Number cv = cp + this->_R;
        //libMesh::Number gamma = cp/cv;
        //libMesh::Number gamma_ratio = gamma/(gamma-1.0);

        libMesh::Number p0 = context.interior_value( this->_p0_var->p0(), qp );

        for (unsigned int i = 0; i != n_p0_dofs; ++i)
          {
            F_p0(i) += (p0/T - this->_p0/this->_T0)*JxW[qp];
            //F_p0(i) -= p0*gamma_ratio*divU*JxW[qp];
          } // End DoF loop i
      }

    return;
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokes<Mu,SH,TC>::assemble_thermo_press_side_time_deriv( bool /*compute_jacobian*/,
                                                                             AssemblyContext& context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_p0_dofs = context.get_dof_indices(this->_p0_var->p0()).size();

    // Element Jacobian * quadrature weight for side integration.
    //const std::vector<libMesh::Real> &JxW_side = context.get_side_fe(this->_temp_vars.T())->get_JxW();

    //const std::vector<Point> &normals = context.get_side_fe(this->_temp_vars.T())->get_normals();

    //libMesh::DenseSubVector<libMesh::Number> &F_p0 = context.get_elem_residual(this->_p0_var->p0()); // residual

    // Physical location of the quadrature points
    //const std::vector<libMesh::Point>& qpoint = context.get_side_fe(this->_temp_vars.T())->get_xyz();

    unsigned int n_qpoints = context.get_side_qrule().n_points();
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        /*
          libMesh::Number T = context.side_value( this->_temp_vars.T(), qp );
          libMesh::Gradient U = ( context.side_value( this->_flow_vars.u(), qp ),
          context.side_value( this->_flow_vars.v(), qp ) );
          libMesh::Gradient grad_T = context.side_gradient( this->_temp_vars.T(), qp );


          libMesh::Number p0 = context.side_value( this->_p0_var->p0(), qp );

          libMesh::Number k = this->_k(T);
          libMesh::Number cp = this->_cp(T);

          libMesh::Number cv = cp + this->_R;
          libMesh::Number gamma = cp/cv;
          libMesh::Number gamma_ratio = gamma/(gamma-1.0);
        */

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
  void LowMachNavierStokes<Mu,SH,TC>::assemble_thermo_press_mass_residual( bool /*compute_jacobian*/,
                                                                           AssemblyContext& context )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_p0_dofs = context.get_dof_indices(this->_p0_var->p0()).size();
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();
    const unsigned int n_p_dofs = context.get_dof_indices(this->_press_var.p()).size();

    // Element Jacobian * quadrature weights for interior integration
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_temp_vars.T())->get_JxW();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_temp_vars.T())->get_phi();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& p_phi =
      context.get_element_fe(this->_press_var.p())->get_phi();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<libMesh::Real> &F_p0 = context.get_elem_residual(this->_p0_var->p0());
    libMesh::DenseSubVector<libMesh::Real> &F_T = context.get_elem_residual(this->_temp_vars.T());
    libMesh::DenseSubVector<libMesh::Real> &F_p = context.get_elem_residual(this->_press_var.p());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
        libMesh::Number T;
        T = context.fixed_interior_value(this->_temp_vars.T(), qp);

        libMesh::Number cp = this->_cp(T);
        libMesh::Number cv = cp + this->_R;
        libMesh::Number gamma = cp/cv;
        libMesh::Number one_over_gamma = 1.0/(gamma-1.0);

        libMesh::Number p0_dot;
        context.interior_rate(this->_p0_var->p0(), qp, p0_dot);

        libMesh::Number p0 = context.fixed_interior_value(this->_p0_var->p0(), qp );

        for (unsigned int i=0; i != n_p0_dofs; i++)
          {
            F_p0(i) -= p0_dot*one_over_gamma*JxW[qp];
          }

        for (unsigned int i=0; i != n_T_dofs; i++)
          {
            F_T(i) += p0_dot*T_phi[i][qp]*JxW[qp];
          }

        for (unsigned int i=0; i != n_p_dofs; i++)
          {
            F_p(i) += p0_dot/p0*p_phi[i][qp]*JxW[qp];
          }

      }
    return;
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokes<Mu,SH,TC>::compute_element_time_derivative_cache( AssemblyContext & context )
  {
    CachedValues & cache = context.get_cached_values();

    const unsigned int n_qpoints = context.get_element_qrule().n_points();

    std::vector<libMesh::Real> u, v, w, T, p, p0;
    u.resize(n_qpoints);
    v.resize(n_qpoints);
    if( this->_flow_vars.dim() > 2 )
      w.resize(n_qpoints);

    T.resize(n_qpoints);
    p.resize(n_qpoints);
    p0.resize(n_qpoints);

    std::vector<libMesh::Gradient> grad_u, grad_v, grad_w, grad_T;
    grad_u.resize(n_qpoints);
    grad_v.resize(n_qpoints);
    if( this->_flow_vars.dim() > 2 )
      grad_w.resize(n_qpoints);

    grad_T.resize(n_qpoints);

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
        u[qp] = context.interior_value(this->_flow_vars.u(), qp);
        v[qp] = context.interior_value(this->_flow_vars.v(), qp);

        grad_u[qp] = context.interior_gradient(this->_flow_vars.u(), qp);
        grad_v[qp] = context.interior_gradient(this->_flow_vars.v(), qp);
        if( this->_flow_vars.dim() > 2 )
          {
            w[qp] = context.interior_value(this->_flow_vars.w(), qp);
            grad_w[qp] = context.interior_gradient(this->_flow_vars.w(), qp);
          }
        T[qp] = context.interior_value(this->_temp_vars.T(), qp);
        grad_T[qp] = context.interior_gradient(this->_temp_vars.T(), qp);

        p[qp] = context.interior_value(this->_press_var.p(), qp);
        p0[qp] = this->get_p0_steady(context, qp);
      }

    cache.set_values(Cache::X_VELOCITY, u);
    cache.set_values(Cache::Y_VELOCITY, v);

    cache.set_gradient_values(Cache::X_VELOCITY_GRAD, grad_u);
    cache.set_gradient_values(Cache::Y_VELOCITY_GRAD, grad_v);

    if(this->_flow_vars.dim() > 2)
      {
        cache.set_values(Cache::Z_VELOCITY, w);
        cache.set_gradient_values(Cache::Z_VELOCITY_GRAD, grad_w);
      }

    cache.set_values(Cache::TEMPERATURE, T);
    cache.set_gradient_values(Cache::TEMPERATURE_GRAD, grad_T);

    cache.set_values(Cache::PRESSURE, p);
    cache.set_values(Cache::THERMO_PRESSURE, p0);
  }

  template<class Mu, class SH, class TC>
  void LowMachNavierStokes<Mu,SH,TC>::compute_postprocessed_quantity( unsigned int quantity_index,
                                                                      const AssemblyContext& context,
                                                                      const libMesh::Point& point,
                                                                      libMesh::Real& value )
  {
    if( quantity_index == this->_rho_index )
      {
        libMesh::Real T = this->T(point,context);
        libMesh::Real p0 = this->get_p0_steady(context,point);

        value = this->rho( T, p0 );
      }

    return;
  }

} // namespace GRINS

// Instantiate
template class GRINS::LowMachNavierStokes<GRINS::ConstantViscosity,GRINS::ConstantSpecificHeat,GRINS::ConstantConductivity>;
