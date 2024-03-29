//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_BOUSSINESQ_BUOYANCY_ADJOINT_STAB_H
#define GRINS_BOUSSINESQ_BUOYANCY_ADJOINT_STAB_H

// GRINS
#include "grins/boussinesq_buoyancy_base.h"
#include "grins/inc_navier_stokes_stab_helper.h"

namespace GRINS
{
  //! Adds Boussinesq bouyancy adjoint stabilization source term
  /*!
    This class implements the adjoint stabilization term for the BoussinesqBuoyancy
    Physics. Intended to be used with IncompressibleNavierStokesAdjointStabilization
    and HeatTransferStabilization.
  */
  template<class Viscosity>
  class BoussinesqBuoyancyAdjointStabilization : public BoussinesqBuoyancyBase
  {
  public:

    BoussinesqBuoyancyAdjointStabilization( const std::string& physics_name, const GetPot& input );

    virtual ~BoussinesqBuoyancyAdjointStabilization() = default;

    virtual void init_context( AssemblyContext & context ) override;

    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext & context ) override;

    virtual void element_constraint( bool compute_jacobian,
                                     AssemblyContext & context ) override;

    // Registers all parameters in this physics and in its property
    // classes
    virtual void register_parameter
    ( const std::string & param_name,
      libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
      const override;

  protected:

    //! Viscosity object
    Viscosity _mu;

    IncompressibleNavierStokesStabilizationHelper _stab_helper;

  };

} // end namespace GRINS
#endif // GRINS_BOUSSINESQ_BUOYANCY_ADJOINT_STAB_H
