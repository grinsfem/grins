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

#ifndef GRINS_INC_NAVIER_STOKES_SPGSM_STAB_H
#define GRINS_INC_NAVIER_STOKES_SPGSM_STAB_H

//GRINS
#include "grins/inc_navier_stokes_stab_base.h"

namespace GRINS
{
  //! Adds VMS-based stabilization to LowMachNavierStokes physics class
  template<class Viscosity>
  class IncompressibleNavierStokesSPGSMStabilization : public IncompressibleNavierStokesStabilizationBase<Viscosity>
  {

  public:

    IncompressibleNavierStokesSPGSMStabilization( const PhysicsName& physics_name, const GetPot& input );

    virtual ~IncompressibleNavierStokesSPGSMStabilization() = default;

    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext& context ) override;

    virtual void element_constraint( bool compute_jacobian,
                                     AssemblyContext & context ) override;

    virtual void mass_residual( bool compute_jacobian,
                                AssemblyContext & context ) override;

  };

} // end namespace GRINS

#endif // GRINS_INC_NAVIER_STOKES_SPGSM_STAB_H
