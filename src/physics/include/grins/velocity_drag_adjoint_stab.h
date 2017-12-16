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


#ifndef GRINS_VELOCITY_DRAG_ADJOINT_STAB_H
#define GRINS_VELOCITY_DRAG_ADJOINT_STAB_H

// GRINS
#include "grins/velocity_drag_base.h"
#include "grins/inc_navier_stokes_stab_helper.h"

// C++
#include <string>

namespace GRINS
{

  //! Physics class for Velocity Drag
  /*
    This physics class imposes a force against the direction of (and
    proportional to an exponent of the magnitude of) a specified
    vector field.
  */
  template<class Viscosity>
  class VelocityDragAdjointStabilization : public VelocityDragBase<Viscosity>
  {
  public:

    VelocityDragAdjointStabilization( const std::string& physics_name, const GetPot& input );

    ~VelocityDragAdjointStabilization();

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

    VelocityDragAdjointStabilization();
  };

} // end namespace block

#endif // GRINS_VELOCITY_DRAG_ADJOINT_STAB_H
