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

// GRINS
#include "grins/reacting_low_mach_navier_stokes_stab_base.h"

#ifndef GRINS_REACTING_LOW_MACH_NAVIER_STOKES_SPGSM_STABILIZATION_H
#define GRINS_REACTING_LOW_MACH_NAVIER_STOKES_SPGSM_STABILIZATION_H

namespace GRINS
{
  template<typename Mixture, typename Evaluator>
  class ReactingLowMachNavierStokesSPGSMStabilization : public ReactingLowMachNavierStokesStabilizationBase<Mixture,Evaluator>
  {
  public:

    ReactingLowMachNavierStokesSPGSMStabilization( const GRINS::PhysicsName& physics_name, const GetPot& input,
                                                   std::unique_ptr<Mixture> & gas_mix);
    virtual ~ReactingLowMachNavierStokesSPGSMStabilization(){};

    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext& context );

    virtual void mass_residual( bool compute_jacobian,
                                AssemblyContext & context );

  private:

    ReactingLowMachNavierStokesSPGSMStabilization();

  };
} // end namespace GRINS

#endif // GRINS_REACTING_LOW_MACH_NAVIER_STOKES_SPGSM_STABILIZATION_H
