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


#ifndef GRINS_SPALART_ALLMARAS_VISCOSITY_H
#define GRINS_SPALART_ALLMARAS_VISCOSITY_H

//GRINS
#include "grins/assembly_context.h"
#include "grins/parameter_user.h"

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/function_base.h"

#include "grins/constant_viscosity.h"
#include "grins/parsed_viscosity.h"
#include "grins/turbulence_fe_variables.h"

#include "libmesh/fem_system.h"

class GetPot;

namespace GRINS
{
  template<class Viscosity>
  class SpalartAllmarasViscosity : public ParameterUser
  {
  public:

    SpalartAllmarasViscosity( const GetPot& input );
    ~SpalartAllmarasViscosity();

    libMesh::Real operator()(AssemblyContext& context, unsigned int qp) const;

    libMesh::Real operator()( const libMesh::Point& p, const libMesh::Real time=0 );

    void init(libMesh::FEMSystem* system);

  protected:

    //! Viscosity object (so we have access to the physical viscosity)
    Viscosity _mu;

    // These are defined for each physics
    TurbulenceFEVariables _turbulence_vars;

  private:

    SpalartAllmarasViscosity();

  };

  /* ------------------------- Inline Functions -------------------------*/
  //inline
  template<class Mu>
  libMesh::Real SpalartAllmarasViscosity<Mu>::operator()(AssemblyContext& context, unsigned int qp) const
  {
    // The physical viscosity
    libMesh::Real mu_physical = this->_mu(context, qp);

    // The unscaled turbulent viscosity (the nu the SA physics solves for)
    libMesh::Real nu = context.interior_value(this->_turbulence_vars.nu_var(),qp);

    // Assert that _mu_value is greater than 0
    if(nu < 0.0)
      {
        libmesh_warning("Negative turbulent viscosity encountered !");

        // We are using a negative S-A model, so will set eddy viscosity to zero
        // if the turbulent viscosity nu < 0.0
        nu = 0.0;
      }

    // Step 1
    libMesh::Real _chi = nu/mu_physical;

    // Step 2
    libMesh::Real _cv1 = 7.1;
    libMesh::Real _fv1 = pow(_chi, 3.0)/(pow(_chi, 3.0) + pow(_cv1, 3.0));

    // Step 3
    libMesh::Real mu_turbulent = nu*_fv1;

    // Compute the value of the total viscosity and return it
    libMesh::Number _mu_value = mu_turbulent + mu_physical; // Turbulent viscosity + physical viscosity

    return _mu_value;
  }

  template<class Mu>
  libMesh::Real SpalartAllmarasViscosity<Mu>::operator()( const libMesh::Point& p, const libMesh::Real time )
  {
    return _mu(p,time);
  }

} // end namespace GRINS

#endif // GRINS_SPALART_ALLMARAS_VISCOSITY_H
