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


#ifndef GRINS_AVERAGED_TURBINE_H
#define GRINS_AVERAGED_TURBINE_H

// GRINS
#include "grins_config.h"
#include "grins/assembly_context.h"
#include "grins/cached_values.h"
#include "grins/inc_navier_stokes_base.h"
#include "grins/averaged_turbine_base.h"

// libMesh
#include "libmesh/getpot.h"

// C++
#include <string>

namespace GRINS
{
  //! Physics class for spatially-averaged turbine
  /*
    This physics class imposes lift/drag forces on velocity as
    affected by a region in which airfoils are moving.  The airfoils
    may also be accelerated or decelerated by external power source or
    sink.
  */
  template<class Viscosity>
  class AveragedTurbine : public AveragedTurbineBase<Viscosity>
  {
  public:

    AveragedTurbine( const std::string& physics_name, const GetPot& input );

    virtual ~AveragedTurbine() = default;

    virtual void init_context( AssemblyContext & context ) override;

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    // Constraint part(s)
    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext & context ) override;

    // Mass residual of the turbine itself
    virtual void nonlocal_mass_residual ( bool compute_jacobian,
                                          AssemblyContext & context ) override;

    // External torque powering the fan or loading the turbine
    virtual void nonlocal_time_derivative ( bool compute_jacobian,
                                            AssemblyContext & context ) override;

  };

} // end namespace block

#endif // GRINS_AVERAGED_TURBINE_H
