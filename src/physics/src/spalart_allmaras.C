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
#include "grins/spalart_allmaras.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/generic_ic_handler.h"
#include "grins/turbulence_models_macro.h"

#include "grins/constant_viscosity.h"
#include "grins/parsed_viscosity.h"
#include "grins/spalart_allmaras_viscosity.h"
#include "grins/variables_parsing.h"
#include "grins/variable_warehouse.h"

// libMesh
#include "libmesh/quadrature.h"
#include "libmesh/elem.h"

namespace GRINS
{

  template<class Mu>
  SpalartAllmaras<Mu>::SpalartAllmaras(const std::string& physics_name, const GetPot& input )
    : TurbulenceModelsBase<Mu>(physics_name, input), // Define class variables
    _flow_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<VelocityVariable>(VariablesParsing::velocity_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
    _press_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PressureFEVariable>(VariablesParsing::press_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
    _turbulence_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<TurbulenceFEVariables>(VariablesParsing::turb_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
    _spalart_allmaras_helper(input),
    _sa_params(input),
    _no_of_walls(input("Physics/"+PhysicsNaming::spalart_allmaras()+"/no_of_walls", 0)),
    _infinite_distance(input("Physics/"+PhysicsNaming::spalart_allmaras()+"/infinite_distance", false))
  {
    _press_var.set_is_constraint_var(true);

    // Loop over the _no_of_walls and fill the wall_ids set
    for(unsigned int i = 0; i != _no_of_walls; i++)
      _wall_ids.insert(input("Physics/"+PhysicsNaming::spalart_allmaras()+"/wall_ids", 0, i ));

    this->_ic_handler = new GenericICHandler( physics_name, input );

    this->check_var_subdomain_consistency(_flow_vars);
    this->check_var_subdomain_consistency(_press_var);
    this->check_var_subdomain_consistency(_turbulence_vars);
  }

  template<class Mu>
  void SpalartAllmaras<Mu>::init_variables( libMesh::FEMSystem* system )
  {
    // Init base class.
    TurbulenceModelsBase<Mu>::init_variables(system);

    // Initialize Boundary Mesh
    this->boundary_mesh.reset(new libMesh::SerialMesh(system->get_mesh().comm() , system->get_mesh().mesh_dimension()) );

    // Use the _wall_ids set to build the boundary mesh object
    (system->get_mesh()).boundary_info->sync(_wall_ids, *boundary_mesh);

    //this->distance_function.reset(new DistanceFunction(system->get_equation_systems(), dynamic_cast<libMesh::UnstructuredMesh&>(system->get_mesh()) ));
    this->distance_function.reset(new DistanceFunction(system->get_equation_systems(), *boundary_mesh));

    // For now, we are hacking this. Without this initialize function being called
    // the distance variable will just be zero. For the channel flow, we are just
    // going to analytically compute the wall distance
    //this->distance_function->initialize();
  }

  template<class Mu>
  void SpalartAllmaras<Mu>::init_context( AssemblyContext& context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    context.get_element_fe(_turbulence_vars.nu())->get_JxW();
    context.get_element_fe(_turbulence_vars.nu())->get_phi();
    context.get_element_fe(_turbulence_vars.nu())->get_dphi();
    context.get_element_fe(_turbulence_vars.nu())->get_xyz();

    context.get_element_fe(_turbulence_vars.nu())->get_phi();
    context.get_element_fe(_turbulence_vars.nu())->get_xyz();

    context.get_side_fe(_turbulence_vars.nu())->get_JxW();
    context.get_side_fe(_turbulence_vars.nu())->get_phi();
    context.get_side_fe(_turbulence_vars.nu())->get_dphi();
    context.get_side_fe(_turbulence_vars.nu())->get_xyz();

    return;
  }

  template<class Mu>
  void SpalartAllmaras<Mu>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Tell the system to march velocity forward in time, but
    // leave p as a constraint only
    system->time_evolving(this->_turbulence_vars.nu(),1);
  }

  template<class Mu>
  void SpalartAllmaras<Mu>::element_time_derivative
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // Get a pointer to the current element, we need this for computing
    // the distance to wall for the  quadrature points
    libMesh::Elem &elem_pointer = context.get_elem();

    // We get some references to cell-specific data that
    // will be used to assemble the linear system.

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_turbulence_vars.nu())->get_JxW();

    // The viscosity shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& nu_phi =
      context.get_element_fe(this->_turbulence_vars.nu())->get_phi();

    // The viscosity shape function gradients (in global coords.)
    // at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& nu_gradphi =
      context.get_element_fe(this->_turbulence_vars.nu())->get_dphi();

    // The number of local degrees of freedom in each variable.
    const unsigned int n_nu_dofs = context.get_dof_indices(this->_turbulence_vars.nu()).size();

    // The subvectors and submatrices we need to fill:
    //
    // K_{\alpha \beta} = R_{\alpha},{\beta} = \partial{ R_{\alpha} } / \partial{ {\beta} } (where R denotes residual)
    // e.g., for \alpha = v and \beta = u we get: K{vu} = R_{v},{u}
    // Note that Kpu, Kpv, Kpw and Fp comes as constraint.

    //libMesh::DenseSubMatrix<libMesh::Number> &Knunu = context.get_elem_jacobian(this->_turbulence_vars.nu(), this->_turbulence_vars.nu()); // R_{nu},{nu}

    libMesh::DenseSubVector<libMesh::Number> &Fnu = context.get_elem_residual(this->_turbulence_vars.nu()); // R_{nu}

    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.get_element_qrule().n_points();

    // Auto pointer to distance fcn evaluated at quad points
    std::unique_ptr< libMesh::DenseVector<libMesh::Real> > distance_qp;

    // Fill the vector of distances to quadrature points
    distance_qp = this->distance_function->interpolate(&elem_pointer, context.get_element_qrule().get_points());

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // Compute the solution & its gradient at the old Newton iterate.
        libMesh::Number nu;
        nu = context.interior_value(this->_turbulence_vars.nu(), qp);

        libMesh::Gradient grad_nu;
        grad_nu = context.interior_gradient(this->_turbulence_vars.nu(), qp);

        libMesh::Real jac = JxW[qp];

        // The physical viscosity
        libMesh::Real mu_qp = this->_mu(context, qp);

        // The vorticity value
        libMesh::Real vorticity_value_qp = this->_spalart_allmaras_helper.vorticity(context, qp);

        // The flow velocity
        libMesh::Number u,v;
        u = context.interior_value(this->_flow_vars.u(), qp);
        v = context.interior_value(this->_flow_vars.v(), qp);

        libMesh::NumberVectorValue U(u,v);
        if (this->_flow_vars.dim() == 3)
          U(2) = context.interior_value(this->_flow_vars.w(), qp);

        //The source term
        libMesh::Real S_tilde = this->_sa_params.source_fn(nu, mu_qp, (*distance_qp)(qp), vorticity_value_qp, _infinite_distance);

        // The ft2 function needed for the negative S-A model
        libMesh::Real chi = nu/mu_qp;
        libMesh::Real f_t2 = this->_sa_params.get_c_t3()*exp(-this->_sa_params.get_c_t4()*chi*chi);

        libMesh::Real source_term = this->_sa_params.get_cb1()*(1 - f_t2)*S_tilde*nu;
        // For a negative turbulent viscosity nu < 0.0 we need to use a different production function
        if(nu < 0.0)
          {
            source_term = this->_sa_params.get_cb1()*(1 - this->_sa_params.get_c_t3())*vorticity_value_qp*nu;
          }

        // The wall destruction term
        libMesh::Real fw = this->_sa_params.destruction_fn(nu, (*distance_qp)(qp), S_tilde, _infinite_distance);

        libMesh::Real nud = 0.0;
        if(_infinite_distance)
          {
            nud = 0.0;
          }
        else
          {
            nud = nu/(*distance_qp)(qp);
          }
        libMesh::Real nud2 = nud*nud;
        libMesh::Real kappa2 = (this->_sa_params.get_kappa())*(this->_sa_params.get_kappa());
        libMesh::Real cw1 = this->_sa_params.get_cb1()/kappa2 + (1.0 + this->_sa_params.get_cb2())/this->_sa_params.get_sigma();
        libMesh::Real destruction_term = (cw1*fw - (this->_sa_params.get_cb1()/kappa2)*f_t2)*nud2;

        // For a negative turbulent viscosity nu < 0.0 we need to use a different production function
        if(nu < 0.0)
          {
            destruction_term = -cw1*nud2;
          }

        libMesh::Real fn1 = 1.0;
        // For a negative turbulent viscosity, fn1 needs to be calculated
        if(nu < 0.0)
          {
            libMesh::Real chi3 = chi*chi*chi;
            fn1 = (this->_sa_params.get_c_n1() + chi3)/(this->_sa_params.get_c_n1() - chi3);
          }

        // First, an i-loop over the viscosity degrees of freedom.
        for (unsigned int i=0; i != n_nu_dofs; i++)
          {
            Fnu(i) += jac *
              ( -this->_rho*(U*grad_nu)*nu_phi[i][qp]  // convection term (assumes incompressibility)
                +source_term*nu_phi[i][qp] // source term
                + (1./this->_sa_params.get_sigma())*(-(mu_qp+(fn1*nu))*grad_nu*nu_gradphi[i][qp] + this->_sa_params.get_cb2()*grad_nu*grad_nu*nu_phi[i][qp]) // diffusion term
                - destruction_term*nu_phi[i][qp]); // destruction term

            // Compute the jacobian if not using numerical jacobians
            if (compute_jacobian)
              {
                libmesh_not_implemented();
              } // end - if (compute_jacobian)

          } // end of the outer dof (i) loop
      } // end of the quadrature point (qp) loop
  }

  template<class K>
  void SpalartAllmaras<K>::mass_residual
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // First we get some references to cell-specific data that
    // will be used to assemble the linear system.

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_turbulence_vars.nu())->get_JxW();

    // The viscosity shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& nu_phi =
      context.get_element_fe(this->_turbulence_vars.nu())->get_phi();

    // The number of local degrees of freedom in each variable.
    const unsigned int n_nu_dofs = context.get_dof_indices(this->_turbulence_vars.nu()).size();

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<libMesh::Real> &F = context.get_elem_residual(this->_turbulence_vars.nu());

    //libMesh::DenseSubMatrix<libMesh::Real> &M = context.get_elem_jacobian(this->_turbulence_vars.nu(), this->_turbulence_vars.nu());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
        // For the mass residual, we need to be a little careful.
        // The time integrator is handling the time-discretization
        // for us so we need to supply M(u_fixed)*u' for the residual.
        // u_fixed will be given by the fixed_interior_value function
        // while u' will be given by the interior_rate function.
        libMesh::Real nu_dot;
        context.interior_rate(this->_turbulence_vars.nu(), qp, nu_dot);

        for (unsigned int i = 0; i != n_nu_dofs; ++i)
          {
            F(i) += -JxW[qp]*this->_rho*nu_dot*nu_phi[i][qp];

            if( compute_jacobian )
              {
                libmesh_not_implemented();
              }// End of check on Jacobian

          } // End of element dof loop

      } // End of the quadrature point loop
  }

  template<class Mu>
  void SpalartAllmaras<Mu>::register_parameter
  ( const std::string & param_name,
    libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
    const
  {
    ParameterUser::register_parameter(param_name, param_pointer);
    this->_mu.register_parameter(param_name, param_pointer);
    this->_sa_params.register_parameter(param_name, param_pointer);
  }

} // namespace GRINS

// Instantiate
INSTANTIATE_TURBULENCE_MODELS_SUBCLASS(SpalartAllmaras);
