//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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

    //! Destructor does not need to delete SharedPtr's.
    virtual ~PhysicsFactory();
    
    //! Builds PhysicsList. This is the primary function of this class.
    GRINS::PhysicsList build(const GetPot& input);

    typedef SharedPtr<Physics> PhysicsPtr;
    typedef std::pair< std::string, PhysicsPtr > PhysicsPair;

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
    /*! The elements in the physics_list may have suffixes attached while
        the required_physics is assumed not to when this function is called.
        So, this helper method will do the check by stripping the suffixes
        from the physics_list element during the search and check if the required
        physics is present.

        \todo This is a linear search, but it really shouldn't matter unless the
              physics_list gets big. */
    bool find_physics( const std::string& required_physics,
                       const GRINS::PhysicsList& physics_list ) const;

    //! Utility function
    void physics_consistency_error( const std::string physics_checked,
				    const std::string physics_required ) const;

  }; // class PhysicsFactory

} // namespace GRINS

#endif // GRINS_PHYSICS_FACTORY_H
