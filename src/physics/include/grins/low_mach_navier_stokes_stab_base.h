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

#ifndef GRINS_LOW_MACH_NAVIER_STOKES_STAB_BASE_H
#define GRINS_LOW_MACH_NAVIER_STOKES_STAB_BASE_H

//GRINS
#include "grins/low_mach_navier_stokes_base.h"
#include "grins/low_mach_navier_stokes_stab_helper.h"

//! GRINS namespace
namespace GRINS
{
  //! Adds VMS-based stabilization to LowMachNavierStokes physics class
  template<class Viscosity, class SpecificHeat, class ThermalConductivity>
  class LowMachNavierStokesStabilizationBase : public LowMachNavierStokesBase<Viscosity,SpecificHeat,ThermalConductivity>
  {

  public:

    LowMachNavierStokesStabilizationBase( const GRINS::PhysicsName& physics_name, const GetPot& input );

    virtual ~LowMachNavierStokesStabilizationBase(){};

    //! Initialize context for added physics variables
    virtual void init_context( AssemblyContext& context );

    libMesh::Real compute_res_continuity_steady( AssemblyContext& context,
                                                 unsigned int qp ) const;

    libMesh::Real compute_res_continuity_transient( AssemblyContext& context,
                                                    unsigned int qp ) const;

    libMesh::RealGradient compute_res_momentum_steady( AssemblyContext& context,
                                                       unsigned int qp ) const;

    libMesh::RealGradient compute_res_momentum_transient( AssemblyContext& context,
                                                          unsigned int qp ) const;

    libMesh::Real compute_res_energy_steady( AssemblyContext& context,
                                             unsigned int qp ) const;

    libMesh::Real compute_res_energy_transient( AssemblyContext& context,
                                                unsigned int qp ) const;

  protected:

    LowMachNavierStokesStabilizationHelper _stab_helper;

  private:

    LowMachNavierStokesStabilizationBase();

  };

} // end namespace GRINS

#endif // GRINS_LOW_MACH_NAVIER_STOKES_STAB_BASE_H
