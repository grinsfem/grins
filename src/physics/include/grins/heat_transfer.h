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

#ifndef GRINS_HEAT_TRANSFER_H
#define GRINS_HEAT_TRANSFER_H

//GRINS
#include "grins/heat_transfer_base.h"

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

  private:

    HeatTransfer();

  };

} // end namespace block

#endif // GRINS_HEAT_TRANSFER_H
