//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - a low Mach number Navier-Stokes Finite-Element Solver
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
#ifndef PHYSICS_SYS_H
#define PHYSICS_SYS_H

#include "config.h"

#include <string>

#include "getpot.h"
#include "libmesh.h"
#include "fem_system.h"
#include "fem_context.h"

// Including this here so every physics class automatically gets it.
#include "variable_name_defaults.h"
#include "var_typedefs.h"

// Boundary conditions
#include "bc_types.h"
#include "boundary_conditions.h"
#include "point_func_base.h"

#ifdef HAVE_GRVY
#include "grvy.h" // GRVY timers
#endif

//! GRINS namespace
namespace GRINS
{
  //! Physics abstract base class. Defines API for physics to be added to MultiphysicsSystem.
  /*!
    This abstract base class defines the API for use within the MultiphysicsSystem. Each physics
    instantiation must initialize its varaibles (and element orders and types), initialize the DiffContext
    for the added variables, and provide element and boundary routines for assembly
    (to be called by MultiphysicsSystem).

    MultiphysicsSystem (through FEMSystem) solves the following equation:

    \f$M(u)\dot{u} = F(u)\f$

    M = mass matrix
    u = solution vector
    F = time derivative

    Note that for the nonlinear system that is solved for implicit
    time stepping is:

    \f$M(u_{\theta})(u^n - u^{n+1}) + \Delta t F(u) = 0\f$
  */

  //TODO: is it F(u) or F(u_{\theta})?

  //  element* routines work on element interiors
  //  side* routines work on element sides

  //  *_time_derivative correspond to calculating terms for F(u)
  //  *_mass_residual correspond to calculating terms for M(u)\dot{u}

  class Physics
  {

  public:

    Physics();
    virtual ~Physics();

    //! Read options from GetPot input file. By default, nothing is read.
    virtual void read_input_options( GetPot& input );

    //! Initialize variables for this physics.
    virtual void init_variables( libMesh::FEMSystem* system ) = 0;

    //! Registers variables for coupled physics.
    /*!
      Each physics might need access to other physics variables, 
      so this method registers them in individual physics. By default,
      no variables are registered.
    */
    virtual void register_variable_indices( GRINS::VariableMap& global_map );

    //! Set which variables are time evolving.
    /*!
      Set those variables which evolve in time (as opposed to variables that behave like constraints).
      This is done separately from init_variables() because the MultiphysicsSystem must initialize
      its base class before time_evolving variables can be set. Default implementation is no
      time evolving variables.
     */
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    //! Initialize context for added physics variables
    virtual void init_context( libMesh::DiffContext &context ) = 0;

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    //! Time dependent part(s) of physics for element interiors
    virtual bool element_time_derivative( bool request_jacobian,
					  libMesh::DiffContext& context,
					  libMesh::FEMSystem* system ) = 0;

    //! Time dependent part(s) of physics for boundaries of elements on the domain boundary
    virtual bool side_time_derivative( bool request_jacobian,
				       libMesh::DiffContext& context,
				       libMesh::FEMSystem* system ) = 0;

    //! Constraint part(s) of physics for element interiors
    virtual bool element_constraint( bool request_jacobian,
				     libMesh::DiffContext& context,
				     libMesh::FEMSystem* system ) = 0;

    //! Constraint part(s) of physics for boundaries of elements on the domain boundary
    virtual bool side_constraint( bool request_jacobian,
				  libMesh::DiffContext& context,
				  libMesh::FEMSystem* system ) = 0;

    //! Mass matrix part(s) for element interiors. All boundary terms lie within the time_derivative part
    virtual bool mass_residual( bool request_jacobian,
				libMesh::DiffContext& context,
				libMesh::FEMSystem* system ) = 0;

    //! Variable map for this physics
    /*
      Force the derived class to build the local physics variable map.
      User must also set _local_variable_map_built = true or the getter
      function will error.
     */
    virtual void build_local_variable_map() = 0;

    //! Returns indices for physics variables
    /*!
      Other physics might need access to this physics variables, so this method returns a copy
      of the map from the std::string name of the variable to the variable index in the system.
    */
    GRINS::VariableMap get_variable_indices_map();

#ifdef USE_GRVY_TIMERS
    void attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer );
#endif

  protected:

    //! Map from std::string variable name to variable index value
    GRINS::VariableMap _var_map;

    //! Need to a way to double check that the local variable map was built
    bool _local_variable_map_built;

    //! Map between boundary id and boundary condition type
    std::map< unsigned int, GRINS::BC_TYPES> _bc_map;

    //! Object that stashes generic boundary condition types
    /** \todo Move this so that only one object is needed. 
	      Perhaps make static? */
    GRINS::BoundaryConditions _bound_conds;

    //! Map between boundary id and boundary functions
    std::map< unsigned int, GRINS::BasePointFuncObj* > _bound_funcs;

#ifdef USE_GRVY_TIMERS
    GRVY::GRVY_Timer_Class* _timer;
#endif

  }; // End Physics class declarations

} // End namespace GRINS

#endif //PHYSICS_SYS_H
