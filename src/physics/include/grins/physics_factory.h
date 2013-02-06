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

#ifndef GRINS_PHYSICS_FACTORY_H
#define GRINS_PHYSICS_FACTORY_H

// C++
#include <string>

// GRINS stuff
#include "grins/var_typedefs.h"

// libMesh forward declarations
class GetPot;

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

    void visc_cond_specheat_error( const std::string& physics,
				   const std::string& conductivity,
				   const std::string& viscosity,
				   const std::string& specific_heat ) const;

  }; // class PhysicsFactory

} // namespace GRINS

#endif // GRINS_PHYSICS_FACTORY_H
