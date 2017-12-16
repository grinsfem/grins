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


#ifndef GRINS_VELOCITY_PENALTY_ADJOINT_STAB_H
#define GRINS_VELOCITY_PENALTY_ADJOINT_STAB_H

// GRINS
#include "grins/velocity_penalty_base.h"
#include "grins/inc_navier_stokes_stab_helper.h"

namespace GRINS
{
  //! Adds Velocity penalty adjoint stabilization source term
  /*!
    This class implements the adjoint stabilization term for the VelocityPenalty
    Physics. Intended to be used with
    IncompressibleNavierStokesAdjointStabilization.
  */
  template<class Viscosity>
  class VelocityPenaltyAdjointStabilization : public VelocityPenaltyBase<Viscosity>
  {
  public:

    VelocityPenaltyAdjointStabilization( const std::string& physics_name, const GetPot& input );

    ~VelocityPenaltyAdjointStabilization();

    virtual void init_context( AssemblyContext& context );

    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext& context );

    virtual void element_constraint( bool compute_jacobian,
                                     AssemblyContext & context );

  protected:

    IncompressibleNavierStokesStabilizationHelper _stab_helper;

  private:

    VelocityPenaltyAdjointStabilization();

  };

} // end namespace GRINS
#endif // GRINS_VELOCITY_PENALTY_ADJOINT_STAB_H
