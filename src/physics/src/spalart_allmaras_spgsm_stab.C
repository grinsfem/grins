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
#include "grins/spalart_allmaras_spgsm_stab.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/constant_viscosity.h"
#include "grins/parsed_viscosity.h"
#include "grins/spalart_allmaras_viscosity.h"
#include "grins/turbulence_models_macro.h"

//libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{

  template<class Mu>
  SpalartAllmarasSPGSMStabilization<Mu>::SpalartAllmarasSPGSMStabilization( const std::string& physics_name,
                                                                            const GetPot& input )
    : SpalartAllmarasStabilizationBase<Mu>(physics_name,input)
  {}

  template<class Mu>
  void SpalartAllmarasSPGSMStabilization<Mu>::init_variables( libMesh::FEMSystem* system )
  {
    // Init base class variables for stab_helper and distance function initialization
    SpalartAllmarasStabilizationBase<Mu>::init_variables(system);

    return;
  }

  template<class Mu>
  void SpalartAllmarasSPGSMStabilization<Mu>::register_parameter
  ( const std::string &param_name, libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer)
    const
  {
    // Register base class parameters
    SpalartAllmarasStabilizationBase<Mu>::register_parameter(param_name, param_pointer);
  }

  template<class Mu>
  void SpalartAllmarasSPGSMStabilization<Mu>::element_time_derivative
  ( bool compute_jacobian,
    AssemblyContext & context )
  {
    // Get a pointer to the current element, we need this for computing the distance to wall for the
    // quadrature points
    libMesh::Elem &elem_pointer = context.get_elem();

    // The number of local degrees of freedom in each variable.
    const unsigned int n_nu_dofs = context.get_dof_indices(this->_turbulence_vars.nu()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_turbulence_vars.nu())->get_JxW();

    // The viscosity shape function gradients (in global coords.)
    // at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& nu_gradphi =
      context.get_element_fe(this->_turbulence_vars.nu())->get_dphi();

    // Quadrature point locations
    //const std::vector<libMesh::Point>& nu_qpoint =
    //context.get_element_fe(this->_turbulence_vars.nu())->get_xyz();

    //libMesh::DenseSubMatrix<libMesh::Number> &Knunu = context.get_elem_jacobian(this->_turbulence_vars.nu(), this->_turbulence_vars.nu()); // R_{nu},{nu}

    libMesh::DenseSubVector<libMesh::Number> &Fnu = context.get_elem_residual(this->_turbulence_vars.nu()); // R_{nu}

    libMesh::FEBase* fe = context.get_element_fe(this->_turbulence_vars.nu());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    // Auto pointer to distance fcn evaluated at quad points
    std::unique_ptr< libMesh::DenseVector<libMesh::Real> > distance_qp;

    // Fill the vector of distances to quadrature points
    distance_qp = this->distance_function->interpolate(&elem_pointer, context.get_element_qrule().get_points());

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Gradient grad_nu;
        grad_nu = context.interior_gradient(this->_turbulence_vars.nu(), qp);

        libMesh::Real jac = JxW[qp];

        // The physical viscosity
        libMesh::Real _mu_qp = this->_mu(context, qp);

        // To be fixed
        // For the channel flow we will just set the distance function analytically
        //(*distance_qp)(qp) = std::min(fabs(y),fabs(1 - y));

        // The flow velocity
        libMesh::Number u,v;
        u = context.interior_value(this->_flow_vars.u(), qp);
        v = context.interior_value(this->_flow_vars.v(), qp);

        libMesh::NumberVectorValue U(u,v);
        if (this->_flow_vars.dim() == 3)
          U(2) = context.interior_value(this->_flow_vars.w(), qp);

        // Stabilization terms

        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::Real tau_spalart = this->_stab_helper.compute_tau_spalart( context, qp, g, G, this->_rho, U, _mu_qp, this->_is_steady );

        libMesh::Number RM_spalart = this->_stab_helper.compute_res_spalart_steady( context, qp, this->_rho, _mu_qp, (*distance_qp)(qp), this->_infinite_distance );

        for (unsigned int i=0; i != n_nu_dofs; i++)
          {
            Fnu(i) += jac*( -tau_spalart*RM_spalart*this->_rho*(U*nu_gradphi[i][qp]) );
          }

        if( compute_jacobian )
          {
            libmesh_not_implemented();
          }

      }
  }

  template<class Mu>
  void SpalartAllmarasSPGSMStabilization<Mu>::mass_residual
  ( bool compute_jacobian, AssemblyContext & context )
  {
    // Get a pointer to the current element, we need this for computing the distance to wall for the
    // quadrature points
    libMesh::Elem &elem_pointer = context.get_elem();

    // The number of local degrees of freedom in each variable.
    const unsigned int n_nu_dofs = context.get_dof_indices(this->_turbulence_vars.nu()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_turbulence_vars.nu())->get_JxW();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& nu_gradphi =
      context.get_element_fe(this->_turbulence_vars.nu())->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> &Fnu = context.get_elem_residual(this->_turbulence_vars.nu()); // R_{nu}

    libMesh::FEBase* fe = context.get_element_fe(this->_turbulence_vars.nu());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    // Auto pointer to distance fcn evaluated at quad points
    std::unique_ptr< libMesh::DenseVector<libMesh::Real> > distance_qp;

    // Fill the vector of distances to quadrature points
    distance_qp = this->distance_function->interpolate(&elem_pointer, context.get_element_qrule().get_points());

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::RealGradient g = this->_stab_helper.compute_g( fe, context, qp );
        libMesh::RealTensor G = this->_stab_helper.compute_G( fe, context, qp );

        libMesh::RealGradient U( context.fixed_interior_value( this->_flow_vars.u(), qp ),
                                 context.fixed_interior_value( this->_flow_vars.v(), qp ) );
        // Compute the viscosity at this qp
        libMesh::Real _mu_qp = this->_mu(context, qp);

        if( this->_flow_vars.dim() == 3 )
          {
            U(2) = context.fixed_interior_value( this->_flow_vars.w(), qp );
          }

        libMesh::Real tau_spalart = this->_stab_helper.compute_tau_spalart( context, qp, g, G, this->_rho, U, _mu_qp, this->_is_steady );

        libMesh::Real RM_spalart = this->_stab_helper.compute_res_spalart_transient( context, qp, this->_rho );

        for (unsigned int i=0; i != n_nu_dofs; i++)
          {
            Fnu(i) += -JxW[qp]*tau_spalart*RM_spalart*this->_rho*(U*nu_gradphi[i][qp]);
          }

        if( compute_jacobian )
          {
            libmesh_not_implemented();
          }

      }
  }

} // end namespace GRINS

// Instantiate
INSTANTIATE_TURBULENCE_MODELS_SUBCLASS(SpalartAllmarasSPGSMStabilization);
