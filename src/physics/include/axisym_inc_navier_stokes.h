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

#ifndef AXISYM_INC_NAVIER_STOKES_H
#define AXISYM_INC_NAVIER_STOKES_H

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
#include "config.h"
#include "physics.h"
#include "pressure_pinning.h"
#include "axisym_inc_navier_stokes_bc_handling.h"

namespace GRINS
{

  //! Physics class for Axisymmetric Incompressible Navier-Stokes
  /*!
    This physics class implements the classical Axisymmetric Incompressible Navier-Stokes equations.
   */
  class AxisymmetricIncompressibleNavierStokes : public Physics
  {
  public:

    AxisymmetricIncompressibleNavierStokes( const std::string& physics_name, const GetPot& input );

    ~AxisymmetricIncompressibleNavierStokes();

    //! Read options from GetPot input file.
    virtual void read_input_options( const GetPot& input );

    //! Initialization of Axisymmetric Navier-Stokes variables
    /*!
      Add velocity and pressure variables to system.
      Note there are only two components of velocity in this case: r and z
     */
    virtual void init_variables( libMesh::FEMSystem* system );

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
    
    //! Handles Dirichlet boundary conditions
    /*! Note that for any generic function specifications, 
      any components not specified will be assigned a zero Dirichlet value. */
    virtual bool side_constraint( bool request_jacobian,
				  libMesh::DiffContext& context,
				  libMesh::FEMSystem* system );

    // Mass matrix part(s)
    virtual bool mass_residual( bool request_jacobian,
				libMesh::DiffContext& context,
				libMesh::FEMSystem* system );

  protected:

    //! Physical dimension of problem
    unsigned int _dim;

    // Index for each owned variable in the system
    //! Index for each r velocity field
    VariableIndex _u_r_var;

    //! Index for each z velocity field
    VariableIndex _u_z_var;

    //! Index for each pressure field
    VariableIndex _p_var;

    // Names of each (owned) variable in the system
    //! r velocity name
    std::string _u_r_var_name;
    
    //! z velocity name
    std::string _u_z_var_name;

    //! pressure name
    std::string _p_var_name;

    //! Element type, read from input
    libMeshEnums::FEFamily _FE_family;

    //! Velocity element order, read from input
    libMeshEnums::Order _V_order;

    //! Pressure element order, read from input
    libMeshEnums::Order _P_order;

    //! Material parameters, read from input
    /** \todo Create objects to allow for function specification */
    libMesh::Number _rho, _mu;

    //! Enable pressure pinning
    bool _pin_pressure;

    PressurePinning _p_pinning;

  private:
    AxisymmetricIncompressibleNavierStokes();

  };

} //End namespace block

#endif // AXISYM_INC_NAVIER_STOKES_H
