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

#ifndef BOUSSINESQ_BUOYANCY_H
#define BOUSSINESQ_BUOYANCY_H

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
  //! Adds Boussinesq bouyancy source term
  /*!
    This class implements the Boussinesq approximation for thermal buoyancy.
    Namely:
    \f$ \mathbf{F} = -\rho_0 \beta_T \left( T - T_0 \right) \mathbf{g} \f$
    where
    \f$ \rho_0 = \f$ reference density, 
    \f$ T_0 = \f$ reference temperature,
    \f$ \beta_T = \f$ coefficient of thermal expansion, and
    \f$ \mathbf{g} = \f$ the gravitional vector.
    This source term to the governing flow equations through the
    element_time_derivative routine. This class requires a flow physics enabled
    and the ConvectiveHeatTransfer physics class enabled.
   */
  class BoussinesqBuoyancy : public Physics
  {
  public:
    
    BoussinesqBuoyancy()
      : Physics()
    {};

    ~BoussinesqBuoyancy()
    {};

    //! Read options from GetPot input file.
    virtual void read_input_options( GetPot& input );

    //! Initialization of BoussinesqBuoyancy variables
    /*!
      There are actually no extra variables
     */
    /*! \todo Perhaps put a default method in the base class so we don't have
      to overload with nothing? */
    virtual void init_variables( libMesh::FEMSystem* system );

    //! Register variables needed by BoussinesqBuoyancy
    /*! This will register the temperature and velocity variables from
      the IncompressibleNavierStokes and ConvectiveHeatTransfer classes.*/
    virtual void register_variable_indices(GRINS::VariableMap &global_map);

    // Context initialization
    /*! Doesn't do anything for BoussinesqBuoyancy since there
      are no new variables registered */
    virtual void init_context( libMesh::DiffContext &context );

    //! Source term contribution for BoussinesqBuoyancy
    /*! This is the main part of the class. This will add the source term to
        the IncompressibleNavierStokes class.
     */
    virtual bool element_time_derivative( bool request_jacobian,
					  libMesh::DiffContext& context,
					  libMesh::FEMSystem* system );

    //! No boundary terms for BoussinesqBuoyancy.
    virtual bool side_time_derivative( bool request_jacobian,
				       libMesh::DiffContext& context,
				       libMesh::FEMSystem* system );

    //! No constraint terms for BoussinesqBuoyancy.
    virtual bool element_constraint( bool request_jacobian,
				     libMesh::DiffContext& context,
				     libMesh::FEMSystem* system );

    //! No boundary terms for BoussinesqBuoyancy.
    virtual bool side_constraint( bool request_jacobian,
				  libMesh::DiffContext& context,
				  libMesh::FEMSystem* system );

    //! No mass terms for BoussinesqBuoyancy.
    virtual bool mass_residual( bool request_jacobian,
				libMesh::DiffContext& context,
				libMesh::FEMSystem* system ); 

    //! No new variables, so no local map
    virtual void build_local_variable_map();

  protected:

    //! Physical dimension of problem
    unsigned int _dim;

    //! Indices for each (registered/non-owned) variable;
    /*!
      This depends on pre-defined set of coupling terms.
     */
    RegtdVariableIndex _u_var; /* Index for x-velocity field */
    RegtdVariableIndex _v_var; /* Index for y-velocity field */
    RegtdVariableIndex _w_var; /* Index for z-velocity field */
    RegtdVariableIndex _T_var; /* Index for Temperature field */

    //! Names of each (non-owned) variable in the system
    std::string _u_var_name, _v_var_name, _w_var_name, _T_var_name;

    //! \f$ \rho_0 = \f$ reference density
    libMesh::Number _rho_ref;

    //! \f$ T_0 = \f$ reference temperature 
    libMesh::Number _T_ref;

    //! \f$ \beta_T = $ coefficient of thermal expansion
    libMesh::Number _beta_T;

    //! Gravitational vector
    Point _g;

  }; // class BoussinesqBuoyancy

} // namespace GRINS
#endif //BOUSSINESQ_BUOYANCY_H
