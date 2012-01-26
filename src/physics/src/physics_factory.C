//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "physics_factory.h"

GRINS::PhysicsFactory::PhysicsFactory( const GetPot& input )
  : _num_physics( input.vector_variable_size("Physics/enabled_physics") )
{
  if( _num_physics < 1 )
    {
      std::cerr << "Error: Must enable at least one physics model" << std::endl;
      libmesh_error();
    }

  // Go through and create a physics object for each physics we're enabling
  for( int i = 0; i < _num_physics; i++ )
    {
      _requested_physics.insert( input("Physics/enabled_physics", "NULL", i ) );
    }

  return;
}

GRINS::PhysicsFactory::~PhysicsFactory()
{
  return;
}

GRINS::PhysicsList GRINS::PhysicsFactory::build()
{
  GRINS::PhysicsList physics_list;

  for( std::set<std::string>::const_iterator physics = _requested_physics.begin();
       physics != _requested_physics.end();
       physics++ )
    {
      this->add_physics( physics, physics_list );
    }

  this->check_physics_consistency( physics_list );

  return physics_list;
}

void GRINS::PhysicsFactory::add_physics( const std::string& physics_to_add,
					 GRINS::PhysicsList& physics_list )
{
  if( physics_to_add == _incompressible_navier_stokes )
    {
      this->physics_list[_incompressible_navier_stokes] = 
	new GRINS::IncompressibleNavierStokes;
    }
  else if( physics_to_add == _axisymmetric_incomp_navier_stokes )
    {
      this->physics_list[_axisymmetric_incomp_navier_stokes] = 
	new GRINS::AxisymmetricIncompNavierStokes;
    }
  else if( physics_to_add == _heat_transfer )
    {
      this->physics_list[_heat_transfer] = new GRINS::HeatTransfer;
    }
  else if( physics_to_add == _axisymmetric_heat_transfer )
    {
      this->physics_list[_axisymmetric_heat_transfer] = 
	new GRINS::AxisymmetricHeatTransfer;
    }
  else if( physics_to_add == _boussinesq_buoyancy )
    {
      this->physics_list[_boussinesq_buoyancy] = new GRINS::BoussinesqBuoyancy;
    }
  else if( physics_to_add == _axisymmetric_boussinesq_buoyancy)
    {
      this->physics_list[_axisymmetric_boussinesq_buoyancy] = 
	new GRINS::AxisymmetricBoussinesqBuoyancy;
    }
  else if( physics_to_add == _axisymmetric_mushy_zone_solidification )
    {
      this->physics_list[_axisymmetric_mushy_zone_solidification] =
	new GRINS::AxisymmetricMushyZoneSolidification;
    }
  else
    {
      std::cerr << "Error: Invalid physics name " << physics_to_add << std::endl;
      libmesh_error();
    }

  return;
}

void GRINS::PhysicsFactory::check_physics_consistency( const GRINS::PhysicsList& physics_list )
{
  for( GRINS::PhysicsListIter physics = physics_list.begin();
       physics != physics_list.end();
       physics++ )
    {
      // For Axisymmetric Incompressible Navier-Stokes, make sure we
      // have two-dimensions
      if( physics->first == _axisymmetric_incomp_navier_stokes )
	{
	  if( this->get_mesh().mesh_dimension() != 2 )
	    {
	      std::cerr << "Error: Dimension of mesh must be 2 for axisymmetric problems."
			<< "Dimension = " 
			<< this->get_mesh().mesh_dimension()
			<< std::endl;
	      libmesh_error();
	    }
	}

      // For HeatTransfer, we need IncompressibleNavierStokes
      if( physics->first == _heat_transfer )
	{
	  if( physics_list.find(_incompressible_navier_stokes) == physics_list.end() )
	    {
	      this->physics_consistency_error( _heat_transfer, _incompressible_navier_stokes  );
	    }
	}

      // For AxisymmetricHeatTransfer, we need AxisymmetricIncompNavierStokes
      if( physics->first == _axisymmetric_heat_transfer )
	{
	  if( physics_list.find(_axisymmetric_incomp_navier_stokes) == 
	      physics_list.end() )
	    {
	      this->physics_consistency_error( _axisymmetric_heat_transfer, 
					       _axisymmetric_incomp_navier_stokes  );
	    }
	}
      // For BoussinesqBuoyancy, we need both HeatTransfer and IncompressibleNavierStokes
      if( physics->first == _boussinesq_buoyancy )
	{
	  if( physics_list.find(_incompressible_navier_stokes) == physics_list.end() )
	    {
	      this->physics_consistency_error( _boussinesq_buoyancy, _incompressible_navier_stokes  );
	    }

	  if( physics_list.find(_heat_transfer) == physics_list.end() )
	    {
	      this->physics_consistency_error( _boussinesq_buoyancy, _heat_transfer  );
	    }
	}

      /* For AxisymmetricBoussinesqBuoyancy, we need both AxisymmetricHeatTransfer 
	 and AxisymmetricIncompNavierStokes */
      if( physics->first == _axisymmetric_boussinesq_buoyancy )
	{
	  if( physics_list.find(_axisymmetric_incomp_navier_stokes) == physics_list.end() )
	    {
	      this->physics_consistency_error( _axisymmetric_boussinesq_buoyancy, 
					       _axisymmetric_incomp_navier_stokes );
	    }

	  if( physics_list.find(_axisymmetric_heat_transfer) == physics_list.end() )
	    {
	      this->physics_consistency_error( _axisymmetric_boussinesq_buoyancy, 
					       _axisymmetric_heat_transfer  );
	    }
	}

      /* For AxisymmetricMushyZoneSolidification, we need both AxisymmetricHeatTransfer 
	 and AxisymmetricIncompNavierStokes */
      if( physics->first == _axisymmetric_mushy_zone_solidification )
	{
	  if( physics_list.find(_axisymmetric_incomp_navier_stokes) == physics_list.end() )
	    {
	      this->physics_consistency_error( _axisymmetric_mushy_zone_solidification, 
					       _axisymmetric_incomp_navier_stokes );
	    }

	  if( physics_list.find(_axisymmetric_heat_transfer) == physics_list.end() )
	    {
	      this->physics_consistency_error( _axisymmetric_mushy_zone_solidification, 
					       _axisymmetric_heat_transfer  );
	    }
	}
    }

  return;
}

void GRINS::PhysicsFactory::physics_consistency_error( const std::string physics_checked,
						       const std::string physics_required )
{
  std::cerr << "Error: " << physics_checked << " physics class requires using "
	    << physics_required << " physics." << std::endl
	    << physics_required << " not found in list of requested physics."
	    << std::endl;

  libmesh_error();	   

  return;
}
