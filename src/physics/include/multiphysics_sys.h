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

#ifndef MULTIPHYSICS_SYS_H
#define MULTIPHYSICS_SYS_H

#include "config.h"

#include <string>

#include "libmesh.h"
#include "boundary_info.h"
#include "fe_base.h"
#include "fe_interface.h"
#include "mesh.h"
#include "quadrature.h"
#include "parameters.h"

// DiffSystem framework files
#include "fem_system.h"
#include "fem_context.h"

#ifdef HAVE_GRVY
// GRVY timers
#include "grvy.h"
#endif

#include "physics.h"

namespace GRINS
{
  //! Interface with libMesh for solving Multiphysics problems.
  /*!
    MultiphysicsSystem (through FEMSystem) solves the following equation:

    \f$M(u)\dot{u} = F(u)\f$
    
    M = mass matrix
    u = solution vector
    F = time derivative

    Note that for the nonlinear system that is solved for implicit
    time stepping is:

    \f$M(u_{\theta})(u^n - u^{n+1}) + \Delta t F(u) = 0\f$ 

    *_time_derivative correspond to calculating terms for \f$F(u)\f$
    *_mass_residual correspond to calculating terms for \f$M(u)\dot{u}\f$
   */
  //TODO: is it F(u) or F(u_{\theta})?
  class MultiphysicsSystem : public FEMSystem
  {    
  public:

    //! Constructor. Will be called by libMesh only.
    MultiphysicsSystem( libMesh::EquationSystems& es,
			const std::string& name,
			const unsigned int number );

    //! Destructor. Clean up all physics allocations.
    ~MultiphysicsSystem();
    
    //! PhysicsList gets built by GRINS::PhysicsFactory and attached here.
    void attach_physics_list( PhysicsList physics_list );

    //! Reads input options for this class and all physics that are enabled
    /*!
      This function reads the input options for the MultiphysicsSystem class and then
      enables each of the requested physics in the system. Finally, the input options
      for each of the physics will be read.
     */
    virtual void read_input_options( const GetPot& input );

    //! System initialization. Calls each physics implementation of init_variables()
    virtual void init_data();

    //! Context initialization. Calls each physics implementation of init_context()
    virtual void init_context( libMesh::DiffContext &context );

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    //! Element interior contributions to \f$F(u)\f$ which have time varying components.
    virtual bool element_time_derivative( bool request_jacobian,
					  libMesh::DiffContext& context );

    //! Boundary contributions to \f$F(u)\f$ which have time varying components.
    virtual bool side_time_derivative( bool request_jacobian,
				       libMesh::DiffContext& context );
    
    //! Element interior contributions to \f$F(u)\f$ which do not have time varying components.
    virtual bool element_constraint( bool request_jacobian,
				     libMesh::DiffContext& context );

    //! Boundary contributions to \f$F(u)\f$ which do not have time varying components.
    virtual bool side_constraint( bool request_jacobian,
				  libMesh::DiffContext& context );
    
    //! Contributions to \f$M(u)\dot{u}\f$
    virtual bool mass_residual( bool request_jacobian,
				libMesh::DiffContext& context ); 

    std::tr1::shared_ptr<GRINS::Physics> get_physics( const std::string physics_name );

#ifdef USE_GRVY_TIMERS
    //! Add GRVY Timer object to system for timing physics.
    void attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer );
#endif    

  private:

    //! Container of pointers to GRINS::Physics classes requested at runtime.
    /*! Set using the attach_physics_list method as construction is taken care
        of by GRINS::PhysicsFactory. */
    PhysicsList _physics_list;

#ifdef USE_GRVY_TIMERS
    GRVY::GRVY_Timer_Class* _timer;
#endif

  };

} //End namespace block

#endif // MULTIPHYSICS_SYS_H
