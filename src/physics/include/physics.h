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
#ifndef PHYSICS_SYS_H
#define PHYSICS_SYS_H

#include "config.h"

#include <string>

//libMesh
#include "getpot.h"
#include "libmesh.h"
#include "fem_system.h"
#include "fem_context.h"

//GRINS
#include "variable_name_defaults.h"
#include "var_typedefs.h"
#include "bc_types.h"
#include "boundary_conditions.h"
#include "grins_physics_names.h"

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

    Physics( const std::string& physics_name );
    virtual ~Physics();

    //! Read options from GetPot input file. By default, nothing is read.
    virtual void read_input_options( const GetPot& input );

    //! Initialize variables for this physics.
    virtual void init_variables( libMesh::FEMSystem* system ) = 0;

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

    virtual void read_bc_data( const GetPot& input );

    virtual int string_to_int( const std::string& bc_type_in );

    virtual void init_bc_data( const GRINS::BoundaryID bc_id, 
			       const std::string& bc_id_string, 
			       const int bc_type, 
			       const GetPot& input );

    void attach_dirichlet_bound_func( GRINS::DBCContainer& dirichlet_bcs );

    void attach_neumann_bound_func( GRINS::NBCContainer& neumann_bcs );

#ifdef USE_GRVY_TIMERS
    void attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer );
#endif

  protected:
    
    //! Name of the physics object. Used for reading physics specific inputs.
    /*! We use a reference because the physics names are const global objects
      in GRINS namespace */
    const std::string& _physics_name;

    //! Map between boundary id and Dirichlet boundary condition type
    std::map< GRINS::BoundaryID, GRINS::BCType> _dirichlet_bc_map;

    //! Map between boundary id and Neumann boundary condition type
    std::map< GRINS::BoundaryID, GRINS::BCType> _neumann_bc_map;

    //! Object that stashes generic boundary condition types
    /** \todo Move this so that only one object is needed. 
	      Perhaps make static? */
    GRINS::BoundaryConditions _bound_conds;

    //! Map between boundary id and general Dirichlet boundary functions
    /*! The user may wish to set a different function for each variable in the physics class,
        for example components of velocity for Navier-Stokes. By design, the user cannot set
	more than 1 function per variable. */
    GRINS::DBCContainer _dirichlet_bound_funcs;
    
    //! Map between boundary id and general Neumann boundary functions
    /*! The user may wish to set a different function for each variable in the physics class. 
        By design, the user cannot set more than 1 function per variable.*/
    GRINS::NBCContainer _neumann_bound_funcs;

#ifdef USE_GRVY_TIMERS
    GRVY::GRVY_Timer_Class* _timer;
#endif

  private:
    Physics();

  }; // End Physics class declarations

} // End namespace GRINS

#endif //PHYSICS_SYS_H
