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

#include "grins_config.h"

#include <string>

//libMesh
#include "getpot.h"
#include "libmesh.h"
#include "fem_system.h"
#include "fem_context.h"
#include "dirichlet_boundaries.h"
#include "zero_function.h"
#include "dof_map.h"

//GRINS
#include "variable_name_defaults.h"
#include "var_typedefs.h"
#include "boundary_conditions.h"
#include "grins_physics_names.h"
#include "grins/bc_handling_base.h"

#ifdef GRINS_HAVE_GRVY
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

    Physics( const GRINS::PhysicsName& physics_name, const GetPot& input );
    virtual ~Physics();

    //! Read options from GetPot input file. By default, nothing is read.
    virtual void read_input_options( const GetPot& input );

    //! Initialize variables for this physics.
    virtual void init_variables( libMesh::FEMSystem* system ) = 0;

    //! Find if current physics is active on supplied element
    virtual bool enabled_on_elem( const libMesh::Elem* elem );

    //! Sets whether this physics is to be solved with a steady solver or not
    /*! Since the member variable is static, only needs to be called on a single
        physics. */
    void set_is_steady( bool is_steady );

    //! Returns whether or not this physics is being solved with a steady solver.
    bool is_steady() const;

    //! Set which variables are time evolving.
    /*!
      Set those variables which evolve in time (as opposed to variables that behave like constraints).
      This is done separately from init_variables() because the MultiphysicsSystem must initialize
      its base class before time_evolving variables can be set. Default implementation is no
      time evolving variables.
     */
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    //! Initialize context for added physics variables
    virtual void init_context( libMesh::FEMContext& context );

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    //! Time dependent part(s) of physics for element interiors
    virtual void element_time_derivative( bool compute_jacobian,
					  libMesh::FEMContext& context );

    //! Time dependent part(s) of physics for boundaries of elements on the domain boundary
    virtual void side_time_derivative( bool compute_jacobian,
				       libMesh::FEMContext& context );

    //! Constraint part(s) of physics for element interiors
    virtual void element_constraint( bool compute_jacobian,
				     libMesh::FEMContext& context );

    //! Constraint part(s) of physics for boundaries of elements on the domain boundary
    virtual void side_constraint( bool compute_jacobian,
				  libMesh::FEMContext& context );

    //! Mass matrix part(s) for element interiors. All boundary terms lie within the time_derivative part
    virtual void mass_residual( bool compute_jacobian,
				libMesh::FEMContext& context );

    void init_bcs( libMesh::FEMSystem* system );

    void attach_neumann_bound_func( GRINS::NBCContainer& neumann_bcs );

    void attach_dirichlet_bound_func( const GRINS::DBCContainer& dirichlet_bc );

#ifdef GRINS_USE_GRVY_TIMERS
    void attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer );
#endif

  protected:
    
    //! Name of the physics object. Used for reading physics specific inputs.
    /*! We use a reference because the physics names are const global objects
      in GRINS namespace */
    const PhysicsName& _physics_name;

    GRINS::BCHandlingBase* _bc_handler;

    //! Subdomains on which the current Physics class is enabled
    std::set<libMesh::subdomain_id_type> _enabled_subdomains;
    
    //! Caches whether or not the solver that's being used is steady or not.
    /*! This is need, for example, in flow stabilization as the tau terms change
        depending on whether the solver is steady or unsteady. */
    static bool _is_steady;

#ifdef GRINS_USE_GRVY_TIMERS
    GRVY::GRVY_Timer_Class* _timer;
#endif

  private:
    Physics();

  }; // End Physics class declarations

} // End namespace GRINS

#endif //PHYSICS_SYS_H
