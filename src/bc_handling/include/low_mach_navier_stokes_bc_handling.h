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
#ifndef LOW_MACH_NAVIER_STOKES_BC_HANDLING_H
#define LOW_MACH_NAVIER_STOKES_BC_HANDLING_H

//libMesh
#include "zero_function.h"

//GRINS
#include "bc_handling_base.h"

namespace GRINS
{
  //! Base class for reading and handling boundary conditions for physics classes
  class LowMachNavierStokesBCHandling : public BCHandlingBase
  {
  public:
    
    LowMachNavierStokesBCHandling( std::string& physics_name, const GetPot& input );
    
    virtual ~LowMachNavierStokesBCHandling();

    virtual int string_to_int( const std::string& bc_type_in ) const;

    virtual void init_bc_data( const GRINS::BoundaryID bc_id, 
			       const std::string& bc_id_string, 
			       const int bc_type, 
			       const GetPot& input );

    void user_init_dirichlet_bcs( libMesh::FEMSystem* system, libMesh::DofMap& dof_map,
				  GRINS::BoundaryID bc_id, GRINS::BCType bc_type ) const;

    
  protected:

    std::string _physics_name;

    std::string _u_var_name, _v_var_name, _w_var_name, _T_var_name;

  private:

    LowMachNavierStokesBCHandling();

    

    enum LMNS_BC_TYPES{NO_SLIP=0, PRESCRIBED_VELOCITY, INFLOW, ISOTHERMAL_WALL,
		       ADIABATIC_WALL, PRESCRIBED_HEAT_FLUX, GENERAL_HEAT_FLUX};

  };
}
#endif // LOW_MACH_NAVIER_STOKES_BC_HANDLING_H
