//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef LOW_MACH_NAVIER_STOKES_BRAACK_STAB_H
#define LOW_MACH_NAVIER_STOKES_BRAACK_STAB_H

//libMesh
#include "time_solver.h"

//GRINS
#include "grins/low_mach_navier_stokes_stab_base.h"

//! GRINS namespace
namespace GRINS
{
  //! Adds VMS-based stabilization to LowMachNavierStokes physics class
  template<class Viscosity, class SpecificHeat, class ThermalConductivity>
  class LowMachNavierStokesBraackStabilization : public LowMachNavierStokesStabilizationBase<Viscosity,SpecificHeat,ThermalConductivity>
  {

  public:

    LowMachNavierStokesBraackStabilization( const GRINS::PhysicsName& physics_name, const GetPot& input );
    virtual ~LowMachNavierStokesBraackStabilization();

    virtual void element_time_derivative( bool compute_jacobian,
					  libMesh::FEMContext& context );

    virtual void mass_residual( bool compute_jacobian,
				libMesh::FEMContext& context );

  protected:

    void assemble_continuity_time_deriv( bool compute_jacobian,
					 libMesh::FEMContext& context );

    void assemble_momentum_time_deriv( bool compute_jacobian,
				       libMesh::FEMContext& context );

    void assemble_energy_time_deriv( bool compute_jacobian,
				     libMesh::FEMContext& context );

    void assemble_continuity_mass_residual( bool compute_jacobian,
					    libMesh::FEMContext& context );

    void assemble_momentum_mass_residual( bool compute_jacobian,
					  libMesh::FEMContext& context );

    void assemble_energy_mass_residual( bool compute_jacobian,
					libMesh::FEMContext& context );
    
  private:
    LowMachNavierStokesBraackStabilization();

  }; // End LowMachNavierStokesBraackStabilization class declarations

} // End namespace GRINS

#endif //LOW_MACH_NAVIER_STOKES_VMS_STAB_H
