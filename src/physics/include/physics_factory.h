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

#ifndef PHYSICS_FACTORY_H
#define PHYSICS_FACTORY_H

#include <string>

// libMesh stuff
#include "getpot.h"

// GRINS stuff
#include "var_typedefs.h"
#include "physics.h"
#include "stokes.h"
#include "inc_navier_stokes.h"
#include "inc_navier_stokes_adjoint_stab.h"
#include "axisym_inc_navier_stokes.h"
#include "heat_transfer.h"
#include "heat_transfer_source.h"
#include "heat_transfer_adjoint_stab.h"
#include "axisym_heat_transfer.h"
#include "boussinesq_buoyancy.h"
#include "axisym_boussinesq_buoyancy.h"
#include "low_mach_navier_stokes.h"
#include "low_mach_navier_stokes_braack_stab.h"
#include "low_mach_navier_stokes_spgsm_stab.h"
#include "low_mach_navier_stokes_vms_stab.h"
#include "grins_physics_names.h"
#include "grins/constant_conductivity.h"
#include "grins/constant_specific_heat.h"
#include "grins/constant_viscosity.h"

namespace GRINS
{
  //! Object for constructing list of physics for simulation
  /*! PhysicsFactory will construct the appropriate GRINS::PhysicsList
      to be handed to the GRINS::MultiphysicsSystem class. The list
      of GRINS::Physics objects constructed is set at run time. 
   */
  class PhysicsFactory
  {
  public:
    
    PhysicsFactory();

    //! Destructor does not need to delete std::tr1::shared_ptr's.
    virtual ~PhysicsFactory();
    
    //! Builds PhysicsList. This is the primary function of this class.
    GRINS::PhysicsList build(const GetPot& input);

  protected:

    //! Figures out which GRINS::Physics pointer to create
    /*! This is the primary method to override if the user wants to extend
        the physics capabilities. The strategy is to conditionally add the
	physics you want, then call the parent PhysicsFactory::add_physics
	function.
     */
    virtual void add_physics( const GetPot& input,
			      const std::string& physics_to_add,
			      GRINS::PhysicsList& physics_list );

    //! Make sure the requested GRINS::Physics classes are consistent
    /*! This is the other method to override (in addition to add_physics)
        for extending physics capabilities. The strategy is to check on
	the physics you've added, then call the parent 
	PhysicsFactory::check_physics_consistency.
     */
    virtual void check_physics_consistency( const GRINS::PhysicsList& physics_list );

    //! Utility function
    void physics_consistency_error( const std::string physics_checked,
				    const std::string physics_required ) const;

    void visc_cond_specheat_error( const std::string& physics, const std::string& conductivity,
				   const std::string& viscosity, const std::string& specific_heat ) const;

  }; // class PhysicsFactory

} // namespace GRINS

#endif //PHYSICS_FACTORY_H
