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

#ifndef AXISYM_BOUSSINESQ_BUOYANCY_H
#define AXISYM_BOUSSINESQ_BUOYANCY_H

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
  //! Adds Axisymmetric Boussinesq bouyancy source term
  /*!
    This class implements the Axisymmetric Boussinesq approximation for thermal buoyancy.
    Namely:
    \f$ \mathbf{F} = -\rho_0 \beta_T \left( T - T_0 \right) \mathbf{g} \f$
    where
    \f$ \rho_0 = \f$ reference density, 
    \f$ T_0 = \f$ reference temperature,
    \f$ \beta_T = \f$ coefficient of thermal expansion, and
    \f$ \mathbf{g} = \f$ the gravitional vector.
    This source term is added to the governing flow equations through the
    element_time_derivative routine. This class requires an axisymmetric flow physics enabled
    and the AxisymmetricHeatTransfer physics class enabled.
   */
  class AxisymmetricBoussinesqBuoyancy : public Physics
  {
  public:
    
    AxisymmetricBoussinesqBuoyancy()
      : Physics()
    {};

    ~AxisymmetricBoussinesqBuoyancy()
    {};

    //! Read options from GetPot input file.
    virtual void read_input_options( GetPot& input );

    //! Initialization of AxisymmetricBoussinesqBuoyancy variables
    /*!
      There are actually no extra variables
     */
    /*! \todo Perhaps put a default method in the base class so we don't have
      to overload with nothing? */
    virtual void init_variables( libMesh::FEMSystem* system );

    //! Register variables needed by AxisymmetricBoussinesqBuoyancy
    /*! This will register the temperature and velocity variables from
      the IncompressibleNavierStokes and ConvectiveHeatTransfer classes.*/
    virtual void register_variable_indices(GRINS::VariableMap &global_map);

    // Context initialization
    /*! Doesn't do anything for AxisymmetricBoussinesqBuoyancy since there
      are no new variables registered */
    virtual void init_context( libMesh::DiffContext &context );

    //! Source term contribution for AxisymmetricBoussinesqBuoyancy
    /*! This is the main part of the class. This will add the source term to
        the AxisymmetricIncompNavierStokes class.
     */
    virtual bool element_time_derivative( bool request_jacobian,
					  libMesh::DiffContext& context,
					  libMesh::FEMSystem* system );

    //! No boundary terms for AxisymmetricBoussinesqBuoyancy.
    virtual bool side_time_derivative( bool request_jacobian,
				       libMesh::DiffContext& context,
				       libMesh::FEMSystem* system );

    //! No constraint terms for AxisymmetricBoussinesqBuoyancy.
    virtual bool element_constraint( bool request_jacobian,
				     libMesh::DiffContext& context,
				     libMesh::FEMSystem* system );

    //! No boundary terms for AxisymmetricBoussinesqBuoyancy.
    virtual bool side_constraint( bool request_jacobian,
				  libMesh::DiffContext& context,
				  libMesh::FEMSystem* system );

    //! No mass terms for AxisymmetricBoussinesqBuoyancy.
    virtual bool mass_residual( bool request_jacobian,
				libMesh::DiffContext& context,
				libMesh::FEMSystem* system ); 

    //! No new variables, so no local map
    virtual void build_local_variable_map();

  protected:

    //! Physical dimension of problem
    unsigned int _dim;

    // Indices for each (registered/non-owned) variable;
    //! Index for registered r-velocity field
    RegtdVariableIndex _u_r_var;

    //! Index for registered z-velocity field
    RegtdVariableIndex _u_z_var;

    //! Index for registered temperature field
    RegtdVariableIndex _T_var;

    // Names of each registered variable in the system

    //! Name of registered r-velocity
    std::string _u_r_var_name;

    //! Name of registered z-velocity
    std::string _u_z_var_name;

    //! Name of registered temperature
    std::string _T_var_name;

    //! \f$ \rho_0 = \f$ reference density
    libMesh::Number _rho_ref;

    //! \f$ T_0 = \f$ reference temperature 
    libMesh::Number _T_ref;

    //! \f$ \beta_T = $ coefficient of thermal expansion
    libMesh::Number _beta_T;

    //! Gravitational vector
    /* \todo This should be stashed in a singleton class and brought in from there */
    Point _g;

  }; // class AxisymmetricBoussinesqBuoyancy

} // namespace GRINS
#endif //AXISYM_BOUSSINESQ_BUOYANCY_H
