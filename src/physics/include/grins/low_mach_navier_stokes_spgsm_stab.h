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

#ifndef GRINS_LOW_MACH_NAVIER_STOKES_SPGSM_STAB_H
#define GRINS_LOW_MACH_NAVIER_STOKES_SPGSM_STAB_H

//GRINS
#include "grins/low_mach_navier_stokes_stab_base.h"

namespace GRINS
{
  //! Adds SPGSM-based stabilization to LowMachNavierStokes physics class
  template<class Viscosity, class SpecificHeat, class ThermalConductivity>
  class LowMachNavierStokesSPGSMStabilization : public LowMachNavierStokesStabilizationBase<Viscosity,SpecificHeat,ThermalConductivity>
  {

  public:

    LowMachNavierStokesSPGSMStabilization( const GRINS::PhysicsName& physics_name, const GetPot& input );
    virtual ~LowMachNavierStokesSPGSMStabilization();

    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext& context );

    virtual void mass_residual( bool compute_jacobian,
                                AssemblyContext & context );

  protected:

    void assemble_continuity_time_deriv( bool compute_jacobian,
                                         AssemblyContext& context );

    void assemble_momentum_time_deriv( bool compute_jacobian,
                                       AssemblyContext& context );

    void assemble_energy_time_deriv( bool compute_jacobian,
                                     AssemblyContext& context );

    void assemble_continuity_mass_residual( bool compute_jacobian,
                                            AssemblyContext& context );

    void assemble_momentum_mass_residual( bool compute_jacobian,
                                          AssemblyContext& context );

    void assemble_energy_mass_residual( bool compute_jacobian,
                                        AssemblyContext& context );

  private:

    LowMachNavierStokesSPGSMStabilization();

  };

} // end namespace GRINS

#endif // GRINS_LOW_MACH_NAVIER_STOKES_SPGSM_STAB_H
