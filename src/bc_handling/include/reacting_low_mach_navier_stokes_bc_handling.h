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
#ifndef GRINS_REACTING_LOW_MACH_NAVIER_STOKES_BC_HANDLING_H
#define GRINS_REACTING_LOW_MACH_NAVIER_STOKES_BC_HANDLING_H

// GRINS
#include "low_mach_navier_stokes_bc_handling.h"

namespace GRINS
{
  class ReactingLowMachNavierStokesBCHandling : public LowMachNavierStokesBCHandling
  {
  public:

    ReactingLowMachNavierStokesBCHandling( const std::string& physics_name, const GetPot& input );

    virtual ~ReactingLowMachNavierStokesBCHandling();

    virtual int string_to_int( const std::string& bc_type_in ) const;

    virtual void init_bc_data( const GRINS::BoundaryID bc_id, 
			       const std::string& bc_id_string, 
			       const int bc_type, 
			       const GetPot& input );

    virtual void user_init_dirichlet_bcs( libMesh::FEMSystem* system, libMesh::DofMap& dof_map,
					  GRINS::BoundaryID bc_id, GRINS::BCType bc_type ) const;

    virtual void init_dirichlet_bcs( libMesh::FEMSystem* system ) const;

  protected:

  private:

    ReactingLowMachNavierStokesBCHandling();

    enum LMNS_BC_TYPES{NO_SLIP=0, PRESCRIBED_VELOCITY, INFLOW, ISOTHERMAL_WALL,
		       ADIABATIC_WALL, PRESCRIBED_HEAT_FLUX, GENERAL_HEAT_FLUX};

    // We need a another container to stash dirichlet values for the speccies
    std::map< GRINS::BoundaryID, Real > _species _values;

    // We also need another map container
    std::map< GRINS::BoundaryID, GRINS::BCType> _temp_bc_map;

  };
}

#endif // GRINS_REACTING_LOW_MACH_NAVIER_STOKES_BC_HANDLING_H
