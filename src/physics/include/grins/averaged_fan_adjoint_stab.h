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


#ifndef GRINS_AVERAGED_FAN_ADJOINT_STAB_H
#define GRINS_AVERAGED_FAN_ADJOINT_STAB_H

// GRINS
#include "grins/averaged_fan_base.h"
#include "grins/inc_navier_stokes_stab_helper.h"

// C++
#include <string>

namespace GRINS
{

  //! Physics class for spatially-averaged fan
  /*
    This physics class imposes lift/drag forces on velocity as
    affected by a region in which airfoils are moving.
  */
  template<class Viscosity>
  class AveragedFanAdjointStabilization : public AveragedFanBase<Viscosity>
  {
  public:

    AveragedFanAdjointStabilization( const std::string& physics_name, const GetPot& input );

    ~AveragedFanAdjointStabilization();

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    virtual void init_context( AssemblyContext& context );

    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext& context );

    virtual void element_constraint( bool compute_jacobian,
                                     AssemblyContext & context );

  protected:

    IncompressibleNavierStokesStabilizationHelper _stab_helper;

  private:

    AveragedFanAdjointStabilization();
  };

} // end namespace block

#endif // GRINS_AVERAGED_FAN_ADJOINT_STAB_H
