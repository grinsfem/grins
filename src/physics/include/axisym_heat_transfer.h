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

#ifndef AXISYM_HEAT_TRANSFER_H
#define AXISYM_HEAT_TRANSFER_H

//libMesh
#include "libmesh.h"
#include "boundary_info.h"
#include "fe_base.h"
#include "fe_interface.h"
#include "mesh.h"
#include "quadrature.h"
#include "parameters.h"
#include "string_to_enum.h"
#include "fem_system.h"
#include "fem_context.h"

//GRINS
#include "grins_config.h"
#include "physics.h"
#include "grins/axisym_heat_transfer_bc_handling.h"

// Conductivity Models
#include "constant_conductivity.h"

namespace GRINS
{

  //! Physics class for Axisymmetric Heat Transfer
  /*
    This physics class implements the classical Axisymmetric Heat Transfer (neglecting viscous dissipation)
   */
  template<class Conductivity>
  class AxisymmetricHeatTransfer : public Physics
  {
  public:

    AxisymmetricHeatTransfer( const std::string& physics_name, const GetPot& input );

    ~AxisymmetricHeatTransfer();

    //! Read options from GetPot input file.
    virtual void read_input_options( const GetPot& input );

    //! Initialization  AxisymmetricHeatTransfer variables
    /*!
      Add velocity and pressure variables to system.
     */
    virtual void init_variables( libMesh::FEMSystem* system );

    //! Sets velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    // Context initialization
    virtual void init_context( libMesh::FEMContext& context );

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    // Time dependent part(s)
    virtual void element_time_derivative( bool compute_jacobian,
					  libMesh::FEMContext& context );

    virtual void side_time_derivative( bool compute_jacobian,
				       libMesh::FEMContext& context );

    // Mass matrix part(s)
    virtual void mass_residual( bool compute_jacobian,
				libMesh::FEMContext& context );

  protected:

    //! Physical dimension of problem
    /*! \todo Make this static member of base class? */
    unsigned int _dim;

    // Indices for each variable;
    //! Index for temperature field
    VariableIndex _T_var;

    //! Index for r-velocity field
    VariableIndex _u_r_var;

    //! Index for z-velocity field
    VariableIndex _u_z_var; 

    // Names of each variable in the system
    //! Name for temperature variable
    std::string _T_var_name;

    //! Name of r-velocity
    std::string _u_r_var_name;

    //! Name of z-velocity
    std::string _u_z_var_name;

    //! Element type, read from input
    libMeshEnums::FEFamily _T_FE_family, _V_FE_family;

    //! Temperature element order, read from input
    libMeshEnums::Order _T_order, _V_order;

    //! Material parameters, read from input
    /*! \todo Need to generalize material parameters. Right now they
              are assumed constant */
    /*! \todo Shouldn't this rho be the same as the one in the flow? Need
              to figure out how to have those shared */
    libMesh::Number _rho, _Cp;

    Conductivity _k;

  private:
    AxisymmetricHeatTransfer();

  };

} //End namespace block

#endif // AXISYM_HEAT_TRANSFER_H
