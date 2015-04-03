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


#ifndef GRINS_STOKES_H
#define GRINS_STOKES_H

//GRINS
#include "grins/pressure_pinning.h"
#include "grins/inc_navier_stokes_bc_handling.h"
#include "grins/inc_navier_stokes_base.h"

namespace GRINS
{

  //! Physics class for Stokes
  /*!
    This physics class implements the classical Stokes equations.
   */
  template<class Viscosity>  
  class Stokes : public IncompressibleNavierStokesBase<Viscosity>
  {
  public:

    Stokes(const std::string& physics_name, const GetPot& input);

    ~Stokes();
    
    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    // Time dependent part(s)
    virtual void element_time_derivative( bool compute_jacobian,
					  AssemblyContext& context,
					  CachedValues& cache );

    // Constraint part(s)
    virtual void element_constraint( bool compute_jacobian,
				     AssemblyContext& context,
				     CachedValues& cache );

    // Mass matrix part(s)
    virtual void mass_residual( bool compute_jacobian,
				AssemblyContext& context,
				CachedValues& cache );

  protected:

    PressurePinning _p_pinning;

    //! Enable pressure pinning
    bool _pin_pressure;
    
  private:

    Stokes();

  };

} // end namespace block

#endif // GRINS_STOKES_H
