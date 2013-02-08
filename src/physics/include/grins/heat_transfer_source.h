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

#ifndef GRINS_HEAT_TRANSFER_SOURCE_H
#define GRINS_HEAT_TRANSFER_SOURCE_H

// GRINS
#include "grins/physics.h"

//libMesh
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"

// libMesh forward declarations
class GetPot;

namespace libMesh
{
  class FEMContext;
  class FEMSystem;
}

namespace GRINS
{  
  //! Adds generic, spatially dependent source term to HeatTransfer physics
  /*! This is templated about the source function. Any suitable source fuction can be
      used so long as its constructor takes a GetPot& and provides and operator( libMesh::Point&) and
      grad( libMesh::Point&) methods which return Real and Gradient respectively.*/
  template< class SourceFunction >
  class HeatTransferSource : public Physics
  {
  public:
    
    HeatTransferSource( const std::string& physics_name, const GetPot& input );

    ~HeatTransferSource();

    virtual void init_variables( libMesh::FEMSystem* system );

    //! Source term contribution for HeatTransferSource
    /*! This is the main part of the class. This will add the source term to
        the HeatTransfer class.
     */
    virtual void element_time_derivative( bool compute_jacobian,
					  libMesh::FEMContext& context,
					  CachedValues& cache );

  protected:

    //! Function that computes source term.
    SourceFunction _source;

    std::string _T_var_name;

    libMeshEnums::FEFamily _T_FE_family;

    libMeshEnums::Order _T_order;

    VariableIndex _T_var; /* Index for temperature field */

  private:
    HeatTransferSource();

  }; // class HeatTransferSource

} // namespace GRINS
#endif // GRINS_HEAT_TRANSFER_SOURCE_H
