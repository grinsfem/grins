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

// This class
#include "grins/physics_factory.h"

// GRINS
#include "grins/cea_thermo.h"
#include "grins/cantera_thermo.h"
#include "grins/constant_transport.h"
#include "grins/cantera_transport.h"
#include "grins/cantera_kinetics.h"
#include "grins/grins_kinetics.h"
#include "grins/physics.h"
#include "grins/stokes.h"
#include "grins/inc_navier_stokes.h"
#include "grins/inc_navier_stokes_adjoint_stab.h"
#include "grins/axisym_inc_navier_stokes.h"
#include "grins/heat_transfer.h"
#include "grins/heat_transfer_source.h"
#include "grins/heat_transfer_adjoint_stab.h"
#include "grins/axisym_heat_transfer.h"
#include "grins/boussinesq_buoyancy.h"
#include "grins/axisym_boussinesq_buoyancy.h"
#include "grins/low_mach_navier_stokes.h"
#include "grins/low_mach_navier_stokes_braack_stab.h"
#include "grins/low_mach_navier_stokes_spgsm_stab.h"
#include "grins/low_mach_navier_stokes_vms_stab.h"
#include "grins/grins_physics_names.h"
#include "grins/constant_conductivity.h"
#include "grins/constant_specific_heat.h"
#include "grins/constant_viscosity.h"
#include "grins/reacting_low_mach_navier_stokes.h"
#include "grins/heat_conduction.h"
#include "grins/constant_source_func.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{

  PhysicsFactory::PhysicsFactory()
  {
    return;
  }

  PhysicsFactory::~PhysicsFactory()
  {
    return;
  }

  PhysicsList PhysicsFactory::build( const GetPot& input )
  {
    PhysicsList physics_list;

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

	for( PhysicsListIter it = physics_list.begin();
	     it != physics_list.end();
	     it++ )
	  {
	    std::cout<< it->first << std::endl;
	  }
	std::cout <<  "==========================================================" << std::endl;
      }

    return physics_list;
  }

  void PhysicsFactory::add_physics( const GetPot& input, 
				    const std::string& physics_to_add,
				    PhysicsList& physics_list )
  {
    typedef std::tr1::shared_ptr<Physics> PhysicsPtr;
    typedef std::pair< std::string, PhysicsPtr > PhysicsPair;

    if( physics_to_add == incompressible_navier_stokes )
      {
	physics_list[physics_to_add] = 
	  PhysicsPtr(new IncompressibleNavierStokes(physics_to_add,input) );
      }
    else if( physics_to_add == stokes )
      {
	physics_list[physics_to_add] =
	  PhysicsPtr(new Stokes(physics_to_add,input));
      }
    else if( physics_to_add == incompressible_navier_stokes_adjoint_stab )
      {
	physics_list[physics_to_add] = 
	  PhysicsPtr(new IncompressibleNavierStokesAdjointStabilization(physics_to_add,input) );
      }
    else if( physics_to_add == axisymmetric_incomp_navier_stokes )
      {
	physics_list[physics_to_add] = 
	  PhysicsPtr(new AxisymmetricIncompressibleNavierStokes(physics_to_add,input) );
      }
    else if( physics_to_add == heat_transfer )
      {
	physics_list[physics_to_add] = 
	  PhysicsPtr(new HeatTransfer(physics_to_add,input));
      }
    else if( physics_to_add == heat_transfer_adjoint_stab )
      {
	physics_list[physics_to_add] = 
	  PhysicsPtr(new HeatTransferAdjointStabilization(physics_to_add,input));
      }
    else if( physics_to_add == heat_transfer_source )
      {
	std::string source_function = input( "Physics/"+physics_to_add+"/source_function", "constant" );
	if( source_function == "constant")
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new HeatTransferSource<ConstantSourceFunction>(physics_to_add,input));
	  }
	else
	  {
	    std::cerr << "================================================================" << std::endl
		      << "Invalid combination of models for " << physics_to_add << std::endl
		      << "Source function  = " << source_function << std::endl
		      << "================================================================" << std::endl;
	    libmesh_error();
	  }
      }
    else if( physics_to_add == axisymmetric_heat_transfer )
      {
	std::string conductivity = input( "Physics/"+axisymmetric_heat_transfer+"/conductivity_model", "constant" );
	if(  conductivity == "constant" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new AxisymmetricHeatTransfer<ConstantConductivity>(physics_to_add,input));
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
	  PhysicsPtr(new BoussinesqBuoyancy(physics_to_add,input));
      }
    else if( physics_to_add == axisymmetric_boussinesq_buoyancy)
      {
	physics_list[physics_to_add] = 
	  PhysicsPtr(new AxisymmetricBoussinesqBuoyancy(physics_to_add,input));
      }
    else if( physics_to_add == "HeatConduction" )
      {
	physics_list[physics_to_add] = 
	  PhysicsPtr(new HeatConduction(physics_to_add,input));
      }
    else if(  physics_to_add == low_mach_navier_stokes )
      {
	std::string conductivity  = input( "Physics/"+low_mach_navier_stokes+"/conductivity_model", "constant" );
	std::string viscosity     = input( "Physics/"+low_mach_navier_stokes+"/viscosity_model", "constant" );
	std::string specific_heat = input( "Physics/"+low_mach_navier_stokes+"/specific_heat_model", "constant" );

	if(  conductivity == "constant" && viscosity == "constant" && specific_heat == "constant" )
	  {
	    physics_list[low_mach_navier_stokes] = 
	      PhysicsPtr(new LowMachNavierStokes<ConstantViscosity,ConstantSpecificHeat,ConstantConductivity>(physics_to_add,input));
	  }
	else
	  {
	    this->visc_cond_specheat_error(physics_to_add, conductivity, viscosity, specific_heat);
	  }
      }
    else if(  physics_to_add == low_mach_navier_stokes_spgsm_stab )
      {
	std::string conductivity  = input( "Physics/"+low_mach_navier_stokes+"/conductivity_model", "constant" );
	std::string viscosity     = input( "Physics/"+low_mach_navier_stokes+"/viscosity_model", "constant" );
	std::string specific_heat = input( "Physics/"+low_mach_navier_stokes+"/specific_heat_model", "constant" );

	if(  conductivity == "constant" && viscosity == "constant" && specific_heat == "constant" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new LowMachNavierStokesSPGSMStabilization<ConstantViscosity,ConstantSpecificHeat,ConstantConductivity>(physics_to_add,input));
	  }
	else
	  {
	    this->visc_cond_specheat_error(physics_to_add, conductivity, viscosity, specific_heat);
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
	      PhysicsPtr(new LowMachNavierStokesVMSStabilization<ConstantViscosity,ConstantSpecificHeat,ConstantConductivity>(physics_to_add,input));
	  }
	else
	  {
	    this->visc_cond_specheat_error(physics_to_add, conductivity, viscosity, specific_heat);
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
	      PhysicsPtr(new LowMachNavierStokesBraackStabilization<ConstantViscosity,ConstantSpecificHeat,ConstantConductivity>(physics_to_add,input));
	  }
	else
	  {
	    this->visc_cond_specheat_error(physics_to_add, conductivity, viscosity, specific_heat);
	  }
      }
    else if( physics_to_add == reacting_low_mach_navier_stokes )
      {
	std::string chem_lib = input( "Physics/"+reacting_low_mach_navier_stokes+"/chemistry_library", "cantera" );
	std::string thermo_lib = input( "Physics/"+reacting_low_mach_navier_stokes+"/thermodynamics_library", "cantera" );
	std::string transport_lib = input( "Physics/"+reacting_low_mach_navier_stokes+"/transport_library", "cantera" );

	if( chem_lib == "cantera" && thermo_lib == "cantera" && transport_lib == "cantera" )
	  {
#ifdef GRINS_HAVE_CANTERA

	    physics_list[physics_to_add] = 
	      PhysicsPtr(new GRINS::ReactingLowMachNavierStokes< GRINS::IdealGasMixture< CanteraThermodynamics,CanteraTransport,CanteraKinetics > >(physics_to_add,input));

#else

	    std::cerr << "Error: Cantera not enable. Cannot use Cantera library."
		      << std::endl;
	    libmesh_error();

#endif // GRINS_HAVE_CANTERA
	  }
	else if( chem_lib == "cantera" && thermo_lib == "cantera" && transport_lib == "grins_constant" )
	  {
#ifdef GRINS_HAVE_CANTERA

	    physics_list[physics_to_add] = 
	      PhysicsPtr(new GRINS::ReactingLowMachNavierStokes< GRINS::IdealGasMixture< CanteraThermodynamics,ConstantTransport,CanteraKinetics > >(physics_to_add,input));

#else

	    std::cerr << "Error: Cantera not enable. Cannot use Cantera library."
		      << std::endl;
	    libmesh_error();

#endif // GRINS_HAVE_CANTERA
	  }
	else if( chem_lib == "cantera" && thermo_lib == "grins_cea" && transport_lib == "grins_constant" )
	  {
#ifdef GRINS_HAVE_CANTERA

	    physics_list[physics_to_add] = 
	      PhysicsPtr(new GRINS::ReactingLowMachNavierStokes< GRINS::IdealGasMixture< CEAThermodynamics,ConstantTransport,CanteraKinetics > >(physics_to_add,input));

#else

	    std::cerr << "Error: Cantera not enable. Cannot use Cantera library."
		      << std::endl;
	    libmesh_error();

#endif // GRINS_HAVE_CANTERA
	  }
	else if( chem_lib == "grins" && thermo_lib == "grins_cea" && transport_lib == "grins_constant" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new GRINS::ReactingLowMachNavierStokes< GRINS::IdealGasMixture< GRINS::CEAThermodynamics,GRINS::ConstantTransport,GRINS::Kinetics > >(physics_to_add,input));
	  }
	else
	  {
	    std::cerr << "Error: Invalid combination of chemistry, transport, and thermodynamics libraries" << std::endl
		      << "       for ReactingLowMachNavierStokes physics." << std::endl
		      << "       chemistry library      = " << chem_lib << std::endl
		      << "       thermodynamics library = " << thermo_lib << std::endl
		      << "       transport library = " << transport_lib << std::endl;
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

  void PhysicsFactory::check_physics_consistency( const PhysicsList& physics_list )
  {
    /*! \todo Need to move the internals of the loop to a separate function that we'll make virtual
      and make this function non-virtual */
    for( PhysicsListIter physics = physics_list.begin();
	 physics != physics_list.end();
	 physics++ )
      {
	// For IncompressibleNavierStokes*Stabilization, we'd better have IncompressibleNavierStokes
	if( physics->first == incompressible_navier_stokes_adjoint_stab )
	  {
	    if( physics_list.find(incompressible_navier_stokes) == physics_list.end() )
	      {
		this->physics_consistency_error( physics->first, incompressible_navier_stokes  );
	      }
	  }

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

	/* For LowMachNavierStokes, there should be nothing else loaded, except
	   for stabilization. */
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

	/* For HeatTransferSource, we'd better have HeatTransfer */
	if( physics->first == heat_transfer_source )
	  {
	    if( physics_list.find(heat_transfer) == physics_list.end() )
	      {
		this->physics_consistency_error( physics->first, heat_transfer  );
	      }
	  }

	/* For HeatTransferAdjointStabilization, we'd better have HeatTransfer */
	if( physics->first == heat_transfer_adjoint_stab )
	  {
	    if( physics_list.find(heat_transfer) == physics_list.end() )
	      {
		this->physics_consistency_error( physics->first, heat_transfer  );
	      }
	  }

	/* For ReactingLowMachNavierStokes, there should be nothing else loaded, except
	   for stabilization. */
	if( physics->first == reacting_low_mach_navier_stokes )
	  {
	    if( physics_list.size() > 2 )
	      {
		std::cerr << "=======================================================" << std::endl
			  << "Error: For physics " << reacting_low_mach_navier_stokes << std::endl
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

  void PhysicsFactory::physics_consistency_error( const std::string physics_checked,
						  const std::string physics_required ) const
  {
    std::cerr << "Error: " << physics_checked << " physics class requires using "
	      << physics_required << " physics." << std::endl
	      << physics_required << " not found in list of requested physics."
	      << std::endl;

    libmesh_error();	   

    return;
  }

  void PhysicsFactory::visc_cond_specheat_error( const std::string& physics,
						 const std::string& conductivity,
						 const std::string& viscosity,
						 const std::string& specific_heat ) const
  {
    std::cerr << "================================================================" << std::endl
	      << "Invalid combination of models for " << physics << std::endl
	      << "Conductivity model  = " << conductivity << std::endl
	      << "Viscosity model     = " << viscosity << std::endl
	      << "Specific heat model = " << specific_heat << std::endl
	      << "================================================================" << std::endl;
    libmesh_error();
  }

} // namespace GRINS
