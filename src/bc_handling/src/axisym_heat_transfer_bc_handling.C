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

#include "axisym_heat_transfer_bc_handling.h"

GRINS::AxisymmetricHeatTransferBCHandling::AxisymmetricHeatTransferBCHandling(std::string& physics_name,
									      const GetPot& input)
  : HeatTransferBCHandling(physics_name, input)
{
  return;
}

GRINS::AxisymmetricHeatTransferBCHandling::~AxisymmetricHeatTransferBCHandling()
{
  return;
}

int GRINS::AxisymmetricHeatTransferBCHandling::string_to_int( const std::string& bc_type ) const
{
  int bc_type_out;

  if( bc_type == "axisymmetrc" )
    bc_type_out = AXISYMMETRIC;
  
  else
    {
      bc_type_out = GRINS::HeatTransferBCHandling::string_to_int( bc_type );
    }

  return bc_type_out;
}

void GRINS::AxisymmetricHeatTransferBCHandling::init_bc_data( const GRINS::BoundaryID bc_id, 
						  const std::string& bc_id_string, 
						  const int bc_type, 
						  const GetPot& input )
{
  switch(bc_type)
    {
    case(AXISYMMETRIC):
      {
	this->set_neumann_bc_type( bc_id, bc_type );
      }
      break;

    default:
      {
	GRINS::HeatTransferBCHandling::init_bc_data( bc_id, bc_id_string, bc_type, input );
      }  
    }// End switch(bc_type)

  return;
}
