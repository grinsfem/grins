//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef HEAT_TRANSFER_H
#define HEAT_TRANSFER_H

//libMesh
#include "libmesh/libmesh.h"
#include "libmesh/boundary_info.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/mesh.h"
#include "libmesh/quadrature.h"
#include "libmesh/parameters.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_system.h"
#include "libmesh/fem_context.h"

//GRINS
#include "grins_config.h"
#include "grins/heat_transfer_base.h"
#include "grins/heat_transfer_bc_handling.h"

namespace GRINS
{

  //! Physics class for Heat Transfer
  /*
    This physics class implements the classical Heat Transfer (neglecting viscous dissipation)
   */
  class HeatTransfer : public HeatTransferBase
  {
  public:

    HeatTransfer( const std::string& physics_name, const GetPot& input );

    ~HeatTransfer();

    //! Read options from GetPot input file.
    virtual void read_input_options( const GetPot& input );
    
    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    // Time dependent part(s)
    virtual void element_time_derivative( bool compute_jacobian,
					  libMesh::FEMContext& context,
					  CachedValues& cache );

    virtual void side_time_derivative( bool compute_jacobian,
				       libMesh::FEMContext& context,
				       CachedValues& cache );

    // Mass matrix part(s)
    virtual void mass_residual( bool compute_jacobian,
				libMesh::FEMContext& context,
				CachedValues& cache );

  protected:
    

  private:
    HeatTransfer();

  };

} //End namespace block

#endif // HEAT_TRANSFER_H
