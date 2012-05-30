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

#include "physics_factory.h"


GRINS::PhysicsFactory::PhysicsFactory()
{
  

  return;
}

GRINS::PhysicsFactory::~PhysicsFactory()
{
  return;
}

GRINS::PhysicsList GRINS::PhysicsFactory::build( const GetPot& input )
{
  GRINS::PhysicsList physics_list;

  int num_physics =  input.vector_variable_size("Physics/enabled_physics");

  if( num_physics < 1 )
    {
      std::cerr << "Error: Must enable at least one physics model" << std::endl;
      libmesh_error();
    }
  
  std::set<std::string> requested_physics;

  // Go through and create a physics object for each physics we're enabling
  for( int i = 0; i < num_physics; i++ )
    {
      requested_physics.insert( input("Physics/enabled_physics", "NULL", i ) );
    }

  for( std::set<std::string>::const_iterator physics = requested_physics.begin();
       physics != requested_physics.end();
       physics++ )
    {
      this->add_physics( input, *physics, physics_list );
    }

  this->check_physics_consistency( physics_list );

  if( input( "screen-options/echo_physics", true ) )
    {
      std::cout << "==========================================================" << std::endl
		<< "List of Enabled Physics:" << std::endl;

      for( GRINS::PhysicsListIter it = physics_list.begin();
	   it != physics_list.end();
	   it++ )
	{
	  std::cout<< it->first << std::endl;
	}
      std::cout <<  "==========================================================" << std::endl;
    }

  return physics_list;
}

void GRINS::PhysicsFactory::add_physics( const GetPot& input, 
					 const std::string& physics_to_add,
					 GRINS::PhysicsList& physics_list )
{
  typedef std::tr1::shared_ptr<GRINS::Physics> PhysicsPtr;
  typedef std::pair< std::string, PhysicsPtr > PhysicsPair;

  if( physics_to_add == incompressible_navier_stokes )
    {
      physics_list[physics_to_add] = 
	PhysicsPtr(new GRINS::IncompressibleNavierStokes(physics_to_add,input) );
    }
  else if( physics_to_add == axisymmetric_incomp_navier_stokes )
    {
      physics_list[physics_to_add] = 
	PhysicsPtr(new GRINS::AxisymmetricIncompressibleNavierStokes(physics_to_add,input) );
    }
  else if( physics_to_add == heat_transfer )
    {
      physics_list[physics_to_add] = 
	PhysicsPtr(new GRINS::HeatTransfer(physics_to_add,input));
    }
  else if( physics_to_add == axisymmetric_heat_transfer )
    {
      std::string conductivity = input( "Physics/"+axisymmetric_heat_transfer+"/conductivity_model", "constant" );
      if(  conductivity == "constant" )
	{
	  physics_list[physics_to_add] = 
	    PhysicsPtr(new GRINS::AxisymmetricHeatTransfer<GRINS::ConstantConductivity>(physics_to_add,input));
	}
      else
	{
	  std::cerr << "Invalid conductivity model " << conductivity << std::endl;
	  libmesh_error();
	}
    }
  else if( physics_to_add == boussinesq_buoyancy )
    {
      physics_list[physics_to_add] = 
	PhysicsPtr(new GRINS::BoussinesqBuoyancy(physics_to_add,input));
    }
  else if( physics_to_add == axisymmetric_boussinesq_buoyancy)
    {
      physics_list[physics_to_add] = 
	PhysicsPtr(new GRINS::AxisymmetricBoussinesqBuoyancy(physics_to_add,input));
    }
  else if(  physics_to_add == low_mach_navier_stokes )
    {
      std::string conductivity  = input( "Physics/"+low_mach_navier_stokes+"/conductivity_model", "constant" );
      std::string viscosity     = input( "Physics/"+low_mach_navier_stokes+"/viscosity_model", "constant" );
      std::string specific_heat = input( "Physics/"+low_mach_navier_stokes+"/specific_heat_model", "constant" );

      if(  conductivity == "constant" && viscosity == "constant" && specific_heat == "constant" )
	{
	  physics_list[low_mach_navier_stokes] = 
	    PhysicsPtr(new GRINS::LowMachNavierStokes<GRINS::ConstantViscosity,GRINS::ConstantSpecificHeat,GRINS::ConstantConductivity>(physics_to_add,input));
	}
      else
	{
	  std::cerr << "================================================================" << std::endl
		    << "Invalid combination of models for " << low_mach_navier_stokes << std::endl
		    << "Conductivity model  = " << conductivity << std::endl
		    << "Viscosity model     = " << viscosity << std::endl
		    << "Specific heat model = " << specific_heat << std::endl
		    << "================================================================" << std::endl;
	  libmesh_error();
	}
    }
  else if(  physics_to_add == low_mach_navier_stokes_vms_stab )
    {
      std::string conductivity  = input( "Physics/"+low_mach_navier_stokes+"/conductivity_model", "constant" );
      std::string viscosity     = input( "Physics/"+low_mach_navier_stokes+"/viscosity_model", "constant" );
      std::string specific_heat = input( "Physics/"+low_mach_navier_stokes+"/specific_heat_model", "constant" );

      if(  conductivity == "constant" && viscosity == "constant" && specific_heat == "constant" )
	{
	  physics_list[physics_to_add] = 
	    PhysicsPtr(new GRINS::LowMachNavierStokesVMSStabilization<GRINS::ConstantViscosity,GRINS::ConstantSpecificHeat,GRINS::ConstantConductivity>(physics_to_add,input));
	}
      else
	{
	  std::cerr << "================================================================" << std::endl
		    << "Invalid combination of models for " << low_mach_navier_stokes_vms_stab << std::endl
		    << "Conductivity model  = " << conductivity << std::endl
		    << "Viscosity model     = " << viscosity << std::endl
		    << "Specific heat model = " << specific_heat << std::endl
		    << "================================================================" << std::endl;
	  libmesh_error();
	}
    }
  else if(  physics_to_add == low_mach_navier_stokes_braack_stab )
    {
      std::string conductivity  = input( "Physics/"+low_mach_navier_stokes+"/conductivity_model", "constant" );
      std::string viscosity     = input( "Physics/"+low_mach_navier_stokes+"/viscosity_model", "constant" );
      std::string specific_heat = input( "Physics/"+low_mach_navier_stokes+"/specific_heat_model", "constant" );

      if(  conductivity == "constant" && viscosity == "constant" && specific_heat == "constant" )
	{
	  physics_list[physics_to_add] = 
	    PhysicsPtr(new GRINS::LowMachNavierStokesBraackStabilization<GRINS::ConstantViscosity,GRINS::ConstantSpecificHeat,GRINS::ConstantConductivity>(physics_to_add,input));
	}
      else
	{
	  std::cerr << "================================================================" << std::endl
		    << "Invalid combination of models for " << low_mach_navier_stokes_vms_stab << std::endl
		    << "Conductivity model  = " << conductivity << std::endl
		    << "Viscosity model     = " << viscosity << std::endl
		    << "Specific heat model = " << specific_heat << std::endl
		    << "================================================================" << std::endl;
	  libmesh_error();
	}
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
  /*! \todo Need to move the internals of the loop to a separate function that we'll make virtual
            and make this function non-virtual */
  for( GRINS::PhysicsListIter physics = physics_list.begin();
       physics != physics_list.end();
       physics++ )
    {
      // For HeatTransfer, we need IncompressibleNavierStokes
      if( physics->first == heat_transfer )
	{
	  if( physics_list.find(incompressible_navier_stokes) == physics_list.end() )
	    {
	      this->physics_consistency_error( heat_transfer, incompressible_navier_stokes  );
	    }
	}

      // For AxisymmetricHeatTransfer, we need AxisymmetricIncompNavierStokes
      if( physics->first == axisymmetric_heat_transfer )
	{
	  if( physics_list.find(axisymmetric_incomp_navier_stokes) == 
	      physics_list.end() )
	    {
	      this->physics_consistency_error( axisymmetric_heat_transfer, 
					       axisymmetric_incomp_navier_stokes  );
	    }
	}
      // For BoussinesqBuoyancy, we need both HeatTransfer and IncompressibleNavierStokes
      if( physics->first == boussinesq_buoyancy )
	{
	  if( physics_list.find(incompressible_navier_stokes) == physics_list.end() )
	    {
	      this->physics_consistency_error( boussinesq_buoyancy, incompressible_navier_stokes  );
	    }

	  if( physics_list.find(heat_transfer) == physics_list.end() )
	    {
	      this->physics_consistency_error( boussinesq_buoyancy, heat_transfer  );
	    }
	}

      /* For AxisymmetricBoussinesqBuoyancy, we need both AxisymmetricHeatTransfer 
	 and AxisymmetricIncompNavierStokes */
      if( physics->first == axisymmetric_boussinesq_buoyancy )
	{
	  if( physics_list.find(axisymmetric_incomp_navier_stokes) == physics_list.end() )
	    {
	      this->physics_consistency_error( axisymmetric_boussinesq_buoyancy, 
					       axisymmetric_incomp_navier_stokes );
	    }

	  if( physics_list.find(axisymmetric_heat_transfer) == physics_list.end() )
	    {
	      this->physics_consistency_error( axisymmetric_boussinesq_buoyancy, 
					       axisymmetric_heat_transfer  );
	    }
	}

      /* For LowMachNavierStokes, there should be nothing else loaded. */
      if( physics->first == low_mach_navier_stokes )
	{
	  if( physics_list.size() > 2 )
	    {
	      std::cerr << "=======================================================" << std::endl
			<< "Error: For physics " << low_mach_navier_stokes << std::endl
			<< "only one stabilization physics is allowed. Detected the" << std::endl
			<< "following:" << std::endl;
	      for( GRINS::PhysicsListIter iter = physics_list.begin();
		   iter != physics_list.end();
		   iter++ )
		{
		  std::cerr << physics->first << std::endl;
		}
	      std::cerr << "=======================================================" << std::endl;
	      libmesh_error();
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
