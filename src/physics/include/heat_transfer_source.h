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

#ifndef HEAT_TRANSFER_SOURCE_H
#define HEAT_TRANSFER_SOURCE_H

// GRINS
#include "config.h"
#include "heat_transfer_base.h"
#include "constant_source_func.h"

namespace GRINS
{  
  //! Adds generic, spatially dependent source term to HeatTransfer physics
  /*! This is templated about the source function. Any suitable source fuction can be
      used so long as its constructor takes a GetPot& and provides and operator( libMesh::Point&) and
      grad( libMesh::Point&) methods which return Real and Gradient respectively.*/
  template< class SourceFunction >
  class HeatTransferSource : public HeatTransferBase
  {
  public:
    
    HeatTransferSource( const std::string& physics_name, const GetPot& input );

    ~HeatTransferSource();

    //! Read options from GetPot input file.
    virtual void read_input_options( const GetPot& input );

    //! Source term contribution for HeatTransferSource
    /*! This is the main part of the class. This will add the source term to
        the HeatTransfer class.
     */
    virtual bool element_time_derivative( bool request_jacobian,
					  libMesh::DiffContext& context,
					  libMesh::FEMSystem* system );

    //! No boundary terms for HeatTransferSource.
    virtual bool side_time_derivative( bool request_jacobian,
				       libMesh::DiffContext& context,
				       libMesh::FEMSystem* system );

    //! No constraint terms for HeatTransferSource.
    virtual bool element_constraint( bool request_jacobian,
				     libMesh::DiffContext& context,
				     libMesh::FEMSystem* system );

    //! No boundary terms for HeatTransferSource.
    virtual bool side_constraint( bool request_jacobian,
				  libMesh::DiffContext& context,
				  libMesh::FEMSystem* system );

    //! No mass terms for HeatTransferSource.
    virtual bool mass_residual( bool request_jacobian,
				libMesh::DiffContext& context,
				libMesh::FEMSystem* system );

  protected:

    //! Function that computes source term.
    SourceFunction _source;

  private:
    HeatTransferSource();

  }; // class HeatTransferSource

} // namespace GRINS
#endif //HEAT_TRANSFER_SOURCE_H
