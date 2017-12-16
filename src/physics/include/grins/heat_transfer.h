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
  template<class Conductivity>
  class HeatTransfer : public HeatTransferBase<Conductivity>
  {
  public:

    HeatTransfer( const std::string& physics_name, const GetPot& input );

    ~HeatTransfer(){};

    //! Register postprocessing variables for HeatTransfer
    virtual void register_postprocessing_vars( const GetPot& input,
                                               PostProcessedQuantities<libMesh::Real>& postprocessing );

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    // Time dependent part(s)
    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext& context );

    // Mass matrix part(s)
    virtual void mass_residual( bool compute_jacobian,
                                AssemblyContext & context );

    //! Compute value of postprocessed quantities at libMesh::Point.
    virtual void compute_postprocessed_quantity( unsigned int quantity_index,
                                                 const AssemblyContext& context,
                                                 const libMesh::Point& point,
                                                 libMesh::Real& value );

  private:

    HeatTransfer();

    //! Index from registering this postprocessed quantity
    unsigned int _k_index;

  };

} // end namespace block

#endif // GRINS_HEAT_TRANSFER_H
