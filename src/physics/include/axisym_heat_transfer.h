//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010,2011 The PECOS Development Team
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

#include "config.h"

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

#include "physics.h"

namespace GRINS
{

  //! Physics class for Axisymmetric Heat Transfer
  /*
    This physics class implements the classical Axisymmetric Heat Transfer (neglecting viscous dissipation)
   */
  class AxisymmetricHeatTransfer : public Physics
  {
  public:

    AxisymmetricHeatTransfer()
      : Physics()
    {};

    ~AxisymmetricHeatTransfer()
    {};

    //! Read options from GetPot input file.
    virtual void read_input_options( GetPot& input );

    //! Initialization  AxisymmetricHeatTransfer variables
    /*!
      Add velocity and pressure variables to system.
     */
    virtual void init_variables( libMesh::FEMSystem* system );

    //! Register external physics variables
    virtual void register_variable_indices( VariableMap& global_map );

    //! Sets velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    // Context initialization
    virtual void init_context( libMesh::DiffContext &context );

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    // Time dependent part(s)
    virtual bool element_time_derivative( bool request_jacobian,
					  libMesh::DiffContext& context,
					  libMesh::FEMSystem* system );

    virtual bool side_time_derivative( bool request_jacobian,
				       libMesh::DiffContext& context,
				       libMesh::FEMSystem* system );

    // Constraint part(s)
    virtual bool element_constraint( bool request_jacobian,
				     libMesh::DiffContext& context,
				     libMesh::FEMSystem* system );

    virtual bool side_constraint( bool request_jacobian,
				  libMesh::DiffContext& context,
				  libMesh::FEMSystem* system );

    // Mass matrix part(s)
    virtual bool mass_residual( bool request_jacobian,
				libMesh::DiffContext& context,
				libMesh::FEMSystem* system );

    virtual void build_local_variable_map();

  protected:

    //! Physical dimension of problem
    /*! \todo Make this static member of base class? */
    unsigned int _dim;

    // Indices for each (owned) variable;
    //! Index for temperature field
    VariableIndex _T_var;

    // Names of each (owned) variable in the system
    //! Name for temperature variable
    std::string _T_var_name;

    // Indices for each (registered/non-owned) variable;
    //! Index for registered r-velocity field
    RegtdVariableIndex _u_r_var;

    //! Index for registered z-velocity field
    RegtdVariableIndex _u_z_var; 

    // Names of each (non-owned) variable in the system
    //! Name of registered r-velocity
    std::string _u_r_var_name;

    //! Name of registered z-velocity
    std::string _u_z_var_name;

    //! Element type, read from input
    libMeshEnums::FEFamily _FE_family;

    //! Temperature element order, read from input
    libMeshEnums::Order _T_order;

    //! Material parameters, read from input
    /*! \todo Need to generalize material parameters. Right now they
              are assumed constant */
    /*! \todo Shouldn't this rho be the same as the one in the flow? Need
              to figure out how to have those shared */
    libMesh::Number _rho, _Cp, _k;

    //! Stash prescribed boundary temperature values
    std::map< unsigned int, double > _T_boundary_values;

    //! Stash prescribed boundary heat flux values
    std::map< unsigned int, double > _q_boundary_values;
  };

} //End namespace block

#endif // AXISYM_HEAT_TRANSFER_H
