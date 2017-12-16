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


#ifndef GRINS_BOUSSINESQ_BUOYANCY_H
#define GRINS_BOUSSINESQ_BUOYANCY_H

// GRINS
#include "grins/boussinesq_buoyancy_base.h"

namespace GRINS
{
  //! Adds Boussinesq bouyancy source term
  /*!
    This class implements the Boussinesq approximation for thermal buoyancy.
    Namely:
    \f$ \mathbf{F} = -\rho_0 \beta_T \left( T - T_0 \right) \mathbf{g} \f$
    where
    \f$ \rho_0 = \f$ reference density,
    \f$ T_0 = \f$ reference temperature,
    \f$ \beta_T = \f$ coefficient of thermal expansion, and
    \f$ \mathbf{g} = \f$ the gravitional vector.
    This source term to the governing flow equations through the
    element_time_derivative routine. This class requires a flow physics enabled
    and the ConvectiveHeatTransfer physics class enabled.
  */
  class BoussinesqBuoyancy : public BoussinesqBuoyancyBase
  {
  public:

    BoussinesqBuoyancy( const std::string& physics_name, const GetPot& input );

    ~BoussinesqBuoyancy();

    //! Source term contribution for BoussinesqBuoyancy
    /*! This is the main part of the class. This will add the source term to
      the IncompressibleNavierStokes class.
    */
    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext& context );

  private:

    BoussinesqBuoyancy();

  };

} // end namespace GRINS
#endif // GRINS_BOUSSINESQ_BUOYANCY_H
