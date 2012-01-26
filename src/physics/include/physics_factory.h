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

#ifndef PHYSICS_FACTORY_H
#define PHYSICS_FACTORY_H

#include <string>

// libMesh stuff
#include "getpot.h"

// GRINS stuff
#include "var_typedefs.h"
#include "physics.h"

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
    
    PhysicsFactory( const GetPot& input );

    //! Destructor does not need to delete AutoPtr's.
    virtual ~PhysicsFactory();
    
    //! Builds PhysicsList. This is the primary function of this class.
    /*! Note that the GRINS::PhysicsList uses libMesh AutoPtr's.
        These are based on std::auto_ptr which means ownership changes.
	\todo Look into Boost install of shared_ptr */
    GRINS::PhysicsList build();

  protected:

    //! Figures out which GRINS::Physics pointer to create
    /*! This is the primary method to override if the user wants to extend
        the physics capabilities. The strategy is to conditionally add the
	physics you want, then call the parent PhysicsFactory::add_physics
	function.
     */
    virtual void add_physics( const std::string& physics_to_add,
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
				    const std::string physics_required )
    
    static const std::string _incompressible_navier_stokes = "IncompressibleNavierStokes";
    static const std::string _axisymmetric_incomp_navier_stokes = "AxisymmetricIncompNavierStokes";
    static const std::string _heat_transfer = "HeatTransfer";
    static const std::string _axisymmetric_heat_transfer = "AxisymmetricHeatTransfer";
    static const std::string _boussinesq_buoyancy = "BoussinesqBuoyancy";
    static const std::string _axisymmetric_boussinesq_buoyancy = "AxisymmetricBoussinesqBuoyancy";
    static const std::string _axisymmetric_mushy_zone_solidification = "AxisymmetricMushyZoneSolidification";

    int _num_physics;
    std::set<std::string> _requested_physics;

  }; // class PhysicsFactory

} // namespace GRINS

#endif //PHYSICS_FACTORY_H
