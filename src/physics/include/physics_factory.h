//-----------------------------------------------------------------------bl-
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
    
    PhysicsFactory();

    //! Destructor does not need to delete AutoPtr's.
    virtual ~PhysicsFactory();
    
    //! Builds PhysicsList. This is the primary function of this class.
    /*! Note that the GRINS::PhysicsList uses libMesh AutoPtr's.
        These are based on std::auto_ptr which means ownership changes.
	\todo Look into Boost install of shared_ptr */
    GRINS::PhysicsList build( const GetPot& input );

  protected:

    //! Figures out which GRINS::Physics pointer to create
    virtual void add_physics( const std::string& physics_to_add );

    //! Make sure the requested GRINS::Physics classes are consistent
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

  }; // class PhysicsFactory

} // namespace GRINS

#endif //PHYSICS_FACTORY_H
