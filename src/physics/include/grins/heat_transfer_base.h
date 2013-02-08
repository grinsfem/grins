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

#ifndef GRINS_HEAT_TRANSFER_BASE_H
#define GRINS_HEAT_TRANSFER_BASE_H

//GRINS
#include "grins/physics.h"

//libMesh
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"

namespace GRINS
{

  //! Physics class for Heat Transfer
  /*
    This physics class implements the classical Heat Transfer (neglecting viscous dissipation)
   */
  class HeatTransferBase : public Physics
  {
  public:

    HeatTransferBase( const std::string& physics_name, const GetPot& input );

    ~HeatTransferBase();

    //! Read options from GetPot input file.
    virtual void read_input_options( const GetPot& input );

    //! Initialization Heat Transfer variables
    /*!
      Add velocity and pressure variables to system.
     */
    virtual void init_variables( libMesh::FEMSystem* system );

    //! Sets velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    // Context initialization
    virtual void init_context( libMesh::FEMContext& context );

  protected:

    //! Physical dimension of problem
    /*! \todo Make this static member of base class? */
    unsigned int _dim;

    //! Indices for each variable;
    VariableIndex _T_var; /* Index for temperature field */
    VariableIndex _u_var; /* Index for x-velocity field */
    VariableIndex _v_var; /* Index for y-velocity field */
    VariableIndex _w_var; /* Index for z-velocity field */

    //! Names of each variable in the system
    std::string _T_var_name;
    std::string _u_var_name, _v_var_name, _w_var_name;

    //! Element type, read from input
    libMeshEnums::FEFamily _T_FE_family, _V_FE_family;

    //! Element orders, read from input
    libMeshEnums::Order _T_order, _V_order;

    //! Material parameters, read from input
    /*! \todo Need to generalize material parameters. Right now they
              are assumed constant */
    /*! \todo Shouldn't this rho be the same as the one in the flow? Need
              to figure out how to have those shared */
    libMesh::Number _rho, _Cp, _k;

  private:
    HeatTransferBase();

  };

} //End namespace block

#endif // GRINS_HEAT_TRANSFER_BASE_H
