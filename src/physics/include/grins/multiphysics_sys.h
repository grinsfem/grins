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

#ifndef GRINS_MULTIPHYSICS_SYS_H
#define GRINS_MULTIPHYSICS_SYS_H

// C++
#include <string>

// GRINS
#include "grins_config.h"
#include "grins/physics.h"

// libMesh
#include "libmesh/fem_system.h"

#ifdef GRINS_HAVE_GRVY
// GRVY timers
#include "grvy.h"
#endif

// libMesh forward declartions
class GetPot;

namespace libMesh
{
  class EquationSystems;
  class DiffContext;
}

namespace GRINS
{
  //! Interface with libMesh for solving Multiphysics problems.
  /*!
    MultiphysicsSystem (through libMesh::FEMSystem) solves the following equation:

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
  class MultiphysicsSystem : public libMesh::FEMSystem
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

    //! Query to check if a particular physics has been enabled
    bool has_physics( const std::string physics_name ) const;

    std::tr1::shared_ptr<GRINS::Physics> get_physics( const std::string physics_name );

    std::tr1::shared_ptr<GRINS::Physics> get_physics( const std::string physics_name ) const;

    void compute_element_cache( const libMesh::FEMContext& context,
				const std::vector<libMesh::Point>& points,
				CachedValues& cache ) const;

#ifdef GRINS_USE_GRVY_TIMERS
    //! Add GRVY Timer object to system for timing physics.
    void attach_grvy_timer( GRVY::GRVY_Timer_Class* grvy_timer );
#endif    

  private:

    //! Container of pointers to GRINS::Physics classes requested at runtime.
    /*! Set using the attach_physics_list method as construction is taken care
        of by GRINS::PhysicsFactory. */
    PhysicsList _physics_list;

    bool _use_numerical_jacobians_only;
    
#ifdef GRINS_USE_GRVY_TIMERS
    GRVY::GRVY_Timer_Class* _timer;
#endif

  };

  inline
  std::tr1::shared_ptr<GRINS::Physics> MultiphysicsSystem::get_physics( const std::string physics_name ) const
  {
    libmesh_assert(_physics_list.find( physics_name ) != _physics_list.end());

    return _physics_list.find(physics_name)->second;
  }

} //End namespace block

#endif // GRINS_MULTIPHYSICS_SYS_H
