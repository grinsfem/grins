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
#ifndef AXISYM_HEAT_TRANSFER_BC_HANDLING_H
#define AXISYM_HEAT_TRANSFER_BC_HANDLING_H

//GRINS
#include "heat_transfer_bc_handling.h"

namespace GRINS
{
  //! Base class for reading and handling boundary conditions for physics classes
  class AxisymmetricHeatTransferBCHandling : public BCHandlingBase
  {
  public:
    
    AxisymmetricHeatTransferBCHandling( const std::string& physics_name, const GetPot& input );
    
    ~AxisymmetricHeatTransferBCHandling();

    int string_to_int( const std::string& bc_type_in ) const;

    void init_bc_data( const GRINS::BoundaryID bc_id, 
		       const std::string& bc_id_string, 
		       const int bc_type, 
		       const GetPot& input );

    void user_apply_neumann_bcs( libMesh::FEMContext& context,
				 GRINS::VariableIndex var,
				 bool request_jacobian,
				 GRINS::BoundaryID bc_id,
				 GRINS::BCType bc_type ) const;
    
    void user_init_dirichlet_bcs( libMesh::FEMSystem* system, 
				  libMesh::DofMap& dof_map,
				  GRINS::BoundaryID bc_id, 
				  GRINS::BCType bc_type ) const;

  private:

    AxisymmetricHeatTransferBCHandling();

    std::string _physics_name;

    std::string _T_var_name;

    enum AHT_BC_TYPES{AXISYMMETRIC=0,
		      ISOTHERMAL_WALL,
		      ADIABATIC_WALL,
		      PRESCRIBED_HEAT_FLUX,
		      GENERAL_HEAT_FLUX};

  };
}
#endif // AXISYM_HEAT_TRANSFER_BC_HANDLING_H
