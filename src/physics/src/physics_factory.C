//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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


// This class
#include "grins/physics_factory.h"

// GRINS
#include "grins/cantera_mixture.h"
#include "grins/cantera_thermo.h"
#include "grins/cantera_transport.h"
#include "grins/cantera_kinetics.h"
#include "grins/cantera_evaluator.h"
#include "grins/physics.h"
#include "grins/stokes.h"
#include "grins/inc_navier_stokes.h"
#include "grins/inc_navier_stokes_adjoint_stab.h"
#include "grins/inc_navier_stokes_spgsm_stab.h"
#include "grins/heat_transfer.h"
#include "grins/heat_transfer_source.h"
#include "grins/heat_transfer_adjoint_stab.h"
#include "grins/heat_transfer_spgsm_stab.h"
#include "grins/axisym_heat_transfer.h"
#include "grins/boussinesq_buoyancy.h"
#include "grins/boussinesq_buoyancy_adjoint_stab.h"
#include "grins/boussinesq_buoyancy_spgsm_stab.h"
#include "grins/axisym_boussinesq_buoyancy.h"
#include "grins/low_mach_navier_stokes.h"
#include "grins/low_mach_navier_stokes_braack_stab.h"
#include "grins/low_mach_navier_stokes_spgsm_stab.h"
#include "grins/low_mach_navier_stokes_vms_stab.h"
#include "grins/averaged_fan.h"
#include "grins/averaged_turbine.h"
#include "grins/scalar_ode.h"
#include "grins/velocity_drag.h"
#include "grins/velocity_penalty.h"
#include "grins/velocity_penalty_adjoint_stab.h"
#include "grins/elastic_membrane.h"
#include "grins/elastic_membrane_constant_pressure.h"
#include "grins/grins_physics_names.h"

#include "grins/spalart_allmaras.h"

#include "grins/constant_conductivity.h"
#include "grins/constant_specific_heat.h"
#include "grins/constant_viscosity.h"
#include "grins/parsed_viscosity.h"
#include "grins/spalart_allmaras_viscosity.h"

#include "grins/reacting_low_mach_navier_stokes.h"
#include "grins/heat_conduction.h"
#include "grins/constant_source_func.h"

#include "grins/antioch_wilke_transport_evaluator.h"
#include "grins/antioch_constant_transport_mixture.h"
#include "grins/antioch_constant_transport_evaluator.h"

#include "grins/hookes_law.h"
#include "grins/incompressible_plane_stress_hyperelasticity.h"
#include "grins/mooney_rivlin.h"

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
    if( physics_to_add == incompressible_navier_stokes )
      {
	std::string viscosity     = input( "Physics/"+incompressible_navier_stokes+"/viscosity_model", "constant" );
	
	if( viscosity == "constant" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new IncompressibleNavierStokes<ConstantViscosity>(physics_to_add,input));
	  }
	else if( viscosity == "parsed" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new IncompressibleNavierStokes<ParsedViscosity>(physics_to_add,input));
	  }
	else if( viscosity == "spalartallmaras" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new IncompressibleNavierStokes<SpalartAllmarasViscosity<ConstantViscosity>>(physics_to_add,input));
	  }
	else
	  {
	    this->visc_error(physics_to_add, viscosity);
	  }          
      }
    else if( physics_to_add == stokes )
      {
	std::string viscosity     = input( "Physics/"+incompressible_navier_stokes+"/viscosity_model", "constant" );
	
	if( viscosity == "constant" )
	  {
	    physics_list[physics_to_add] =
	      PhysicsPtr(new Stokes<ConstantViscosity>(physics_to_add,input));
	  }
	else if( viscosity == "parsed" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new Stokes<ParsedViscosity>(physics_to_add,input));
	  }
	else if( viscosity == "spalartallmaras" ) // For model adaptive cases
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new Stokes<SpalartAllmarasViscosity<ConstantViscosity>>(physics_to_add,input));
	  }	
	else
	  {
	    this->visc_error(physics_to_add, viscosity);
	  }
      }      
    else if( physics_to_add == incompressible_navier_stokes_adjoint_stab )
      {
	std::string viscosity     = input( "Physics/"+incompressible_navier_stokes+"/viscosity_model", "constant" );

	if( viscosity == "constant" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new IncompressibleNavierStokesAdjointStabilization<ConstantViscosity>(physics_to_add,input) );
	  }
	else if( viscosity == "parsed" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new IncompressibleNavierStokesAdjointStabilization<ParsedViscosity>(physics_to_add,input));
	  }
	else if( viscosity == "spalartallmaras" ) 
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new IncompressibleNavierStokesAdjointStabilization<SpalartAllmarasViscosity<ConstantViscosity>>(physics_to_add,input));
	  }		
	else
	  {
	    this->visc_error(physics_to_add, viscosity);
	  }  
      }
    else if( physics_to_add == incompressible_navier_stokes_spgsm_stab )
      {
	std::string viscosity     = input( "Physics/"+incompressible_navier_stokes+"/viscosity_model", "constant" );
	
	if( viscosity == "constant" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new IncompressibleNavierStokesSPGSMStabilization<ConstantViscosity>(physics_to_add,input) );
	  }
	else if( viscosity == "parsed" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new IncompressibleNavierStokesSPGSMStabilization<ParsedViscosity>(physics_to_add,input));
	  }
	else if( viscosity == "spalartallmaras" ) 
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new IncompressibleNavierStokesSPGSMStabilization<SpalartAllmarasViscosity<ConstantViscosity>>(physics_to_add,input));
	  }		
	else
	  {
	    this->visc_error(physics_to_add, viscosity);
	  }  
      }
    else if( physics_to_add == velocity_drag )
      {
	std::string viscosity     = input( "Physics/"+incompressible_navier_stokes+"/viscosity_model", "constant" );

	if( viscosity == "constant" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new VelocityDrag<ConstantViscosity>(physics_to_add,input));
	  }
	else if( viscosity == "parsed" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new VelocityDrag<ParsedViscosity>(physics_to_add,input));
	  }
	else if( viscosity == "spalartallmaras" ) 
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new VelocityDrag<SpalartAllmarasViscosity<ConstantViscosity>>(physics_to_add,input));
	  }			
	else
	  {
	    this->visc_error(physics_to_add, viscosity);
	  }
      }
    else if( physics_to_add == velocity_penalty )
      {
	std::string viscosity     = input( "Physics/"+incompressible_navier_stokes+"/viscosity_model", "constant" );

	if( viscosity == "constant" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new VelocityPenalty<ConstantViscosity>(physics_to_add,input));
	  }
	else if( viscosity == "parsed" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new VelocityPenalty<ParsedViscosity>(physics_to_add,input));
	  }
	else if( viscosity == "spalartallmaras" ) 
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new VelocityPenalty<SpalartAllmarasViscosity<ConstantViscosity>>(physics_to_add,input));
	  }			
	else
	  {
	    this->visc_error(physics_to_add, viscosity);
	  }
      }
    else if( physics_to_add == velocity_penalty_adjoint_stab )
      {
	std::string viscosity     = input( "Physics/"+incompressible_navier_stokes+"/viscosity_model", "constant" );

	if( viscosity == "constant" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new VelocityPenaltyAdjointStabilization<ConstantViscosity>(physics_to_add,input));
	  }
	else if( viscosity == "parsed" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new VelocityPenaltyAdjointStabilization<ParsedViscosity>(physics_to_add,input));
	  }
	else if( viscosity == "spalartallmaras" ) 
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new VelocityPenaltyAdjointStabilization<SpalartAllmarasViscosity<ConstantViscosity>>(physics_to_add,input));
	  }			
	else
	  {
	    this->visc_error(physics_to_add, viscosity);
	  }
      }
    else if( physics_to_add == averaged_fan )
      {
	std::string viscosity     = input( "Physics/"+incompressible_navier_stokes+"/viscosity_model", "constant" );

	if( viscosity == "constant" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new AveragedFan<ConstantViscosity>(physics_to_add,input));
	  }
	else if( viscosity == "parsed" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new AveragedFan<ParsedViscosity>(physics_to_add,input));
	  }
	else if( viscosity == "spalartallmaras" ) 
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new AveragedFan<SpalartAllmarasViscosity<ConstantViscosity>>(physics_to_add,input));
	  }			
	else
	  {
	    this->visc_error(physics_to_add, viscosity);
	  }
      }
    else if( physics_to_add == averaged_turbine )
      {
	std::string viscosity     = input( "Physics/"+incompressible_navier_stokes+"/viscosity_model", "constant" );

	if( viscosity == "constant" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new AveragedTurbine<ConstantViscosity>(physics_to_add,input));
	  }
	else if( viscosity == "parsed" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new AveragedTurbine<ParsedViscosity>(physics_to_add,input));
	  }
	else if( viscosity == "spalartallmaras" ) 
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new AveragedTurbine<SpalartAllmarasViscosity<ConstantViscosity>>(physics_to_add,input));
	  }			
	else
	  {
	    this->visc_error(physics_to_add, viscosity);
	  }
      }
    else if( physics_to_add == spalart_allmaras )
      {
	std::string viscosity     = input( "Physics/"+incompressible_navier_stokes+"/viscosity_model", "constant" );

	if( viscosity == "spalartallmaras" ) 
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new AveragedTurbine<SpalartAllmarasViscosity<ConstantViscosity>>(physics_to_add,input));
	  }			
	else     //Viscosity has to be SA viscosity if SA turbulence model is being used
	  {
	    this->visc_error(physics_to_add, viscosity);
	  }
      }
    else if( physics_to_add == scalar_ode )
      {
	physics_list[physics_to_add] = 
	  PhysicsPtr(new ScalarODE(physics_to_add,input));
      }
    else if( physics_to_add == heat_transfer )
      {
	std::string conductivity     = input( "Physics/"+heat_transfer+"/conductivity_model", "constant" );
	
	if( conductivity == "constant" )
	  {	    
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new HeatTransfer<ConstantConductivity>(physics_to_add,input));
	  }
	else if( conductivity == "parsed" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new HeatTransfer<ParsedConductivity>(physics_to_add,input));
	  }
	else
	  {
	    this->conductivity_error(physics_to_add, conductivity);
	  }
      }
    else if( physics_to_add == heat_transfer_adjoint_stab )
      {
	std::string conductivity     = input( "Physics/"+heat_transfer+"/conductivity_model", "constant" );
	
	if( conductivity == "constant" )
	  {	    
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new HeatTransferAdjointStabilization<ConstantConductivity>(physics_to_add,input));
	  }
	else if( conductivity == "parsed" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new HeatTransferAdjointStabilization<ParsedConductivity>(physics_to_add,input));
	  }
	else
	  {
	    this->conductivity_error(physics_to_add, conductivity);
	  }		
      }
    else if( physics_to_add == heat_transfer_spgsm_stab )
      {
	std::string conductivity     = input( "Physics/"+heat_transfer+"/conductivity_model", "constant" );
	
	if( conductivity == "constant" )
	  {	    
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new HeatTransferSPGSMStabilization<ConstantConductivity>(physics_to_add,input));
	  }
	else if( conductivity == "parsed" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new HeatTransferSPGSMStabilization<ParsedConductivity>(physics_to_add,input));
	  }
	else
	  {
	    this->conductivity_error(physics_to_add, conductivity);
	  }        
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
    else if( physics_to_add == boussinesq_buoyancy_adjoint_stab )
      {
	physics_list[physics_to_add] = 
	  PhysicsPtr(new BoussinesqBuoyancyAdjointStabilization(physics_to_add,input));
      }
    else if( physics_to_add == boussinesq_buoyancy_spgsm_stab )
      {
        physics_list[physics_to_add] =
          PhysicsPtr(new BoussinesqBuoyancySPGSMStabilization(physics_to_add,input));
      }
    else if( physics_to_add == axisymmetric_boussinesq_buoyancy)
      {
	physics_list[physics_to_add] = 
	  PhysicsPtr(new AxisymmetricBoussinesqBuoyancy(physics_to_add,input));
      }
    else if( physics_to_add == "HeatConduction" )
      {
	std::string conductivity     = input( "Physics/"+heat_transfer+"/conductivity_model", "constant" );
	
	if( conductivity == "constant" )
	  {	    
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new HeatConduction<ConstantConductivity>(physics_to_add,input));
	  }
	else if( conductivity == "parsed" )
	  {
	    physics_list[physics_to_add] = 
	      PhysicsPtr(new HeatConduction<ParsedConductivity>(physics_to_add,input));
	  }
	else
	  {
	    this->conductivity_error(physics_to_add, conductivity);
	  }	
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
        this->add_reacting_low_mach( input, physics_to_add, physics_list );
      }
    else if( physics_to_add == elastic_membrane )
      {
        std::string elasticity_model = input("Physics/"+elastic_membrane+"/elasticity_model", "HookesLaw" );

        if( elasticity_model == std::string("HookesLaw") )
          {
            physics_list[physics_to_add] =
              // We need to track \lambda as an indendent variable
              PhysicsPtr(new ElasticMembrane<HookesLaw>(physics_to_add,input,true));
          }
        else if( elasticity_model == std::string("MooneyRivlin") )
          {
            physics_list[physics_to_add] =
              // \lambda determined from incompressiblity
              PhysicsPtr(new ElasticMembrane<IncompressiblePlaneStressHyperelasticity<MooneyRivlin> >(physics_to_add,input,false));
          }
        else
          {
            std::cerr << "Error: Invalid elasticity_model: " << elasticity_model << std::endl
                      << "       Valid selections are: Hookean" << std::endl;
            libmesh_error();
          }
      }
    else if( physics_to_add == elastic_membrane_constant_pressure )
      {
        physics_list[physics_to_add] =
          PhysicsPtr(new ElasticMembraneConstantPressure(physics_to_add,input));
      }
    else
      {
        std::cerr << "Error: Invalid physics name " << physics_to_add << std::endl;
        libmesh_error();
      }
  
    return;
  }

  void PhysicsFactory::add_reacting_low_mach( const GetPot& input,
                                              const std::string& physics_to_add,
                                              GRINS::PhysicsList& physics_list )
  {
    std::string thermochem_lib = input( "Physics/"+reacting_low_mach_navier_stokes+"/thermochemistry_library", "DIE!" );;

    if( thermochem_lib == "cantera" )
      {
#ifdef GRINS_HAVE_CANTERA
        physics_list[physics_to_add] = 
          PhysicsPtr(new GRINS::ReactingLowMachNavierStokes<CanteraMixture,CanteraEvaluator>(physics_to_add,input));
#else
        std::cerr << "Error: Cantera not enabled. Cannot use Cantera library."
                  << std::endl;
        libmesh_error();

#endif // GRINS_HAVE_CANTERA
      }
    else if( thermochem_lib == "antioch" )
      {
#ifdef GRINS_HAVE_ANTIOCH
        std::string mixing_model = input( "Physics/Antioch/mixing_model" , "wilke" );

        std::string thermo_model = input( "Physics/Antioch/thermo_model", "stat_mech");
        std::string viscosity_model = input( "Physics/Antioch/viscosity_model", "blottner");
        std::string conductivity_model = input( "Physics/Antioch/conductivity_model", "eucken");
        std::string diffusivity_model = input( "Physics/Antioch/diffusivity_model", "constant_lewis");

        if( mixing_model == std::string("wilke") )
          {
            if( (thermo_model == std::string("stat_mech")) &&
                (diffusivity_model == std::string("constant_lewis")) &&
                (conductivity_model == std::string("eucken")) &&
                (viscosity_model == std::string("sutherland")) )
              {
                physics_list[physics_to_add] = 
                  PhysicsPtr(new GRINS::ReactingLowMachNavierStokes<
                             GRINS::AntiochWilkeTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                 Antioch::MixtureViscosity<Antioch::SutherlandViscosity<libMesh::Real> >,
                                                                 Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                 Antioch::ConstantLewisDiffusivity<libMesh::Real> >,
                             GRINS::AntiochWilkeTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                  Antioch::MixtureViscosity<Antioch::SutherlandViscosity<libMesh::Real> >,
                                                                  Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                  Antioch::ConstantLewisDiffusivity<libMesh::Real> > >(physics_to_add,input) );
              }
            else if( (thermo_model == std::string("stat_mech")) &&
                     (diffusivity_model == std::string("constant_lewis")) &&
                     (conductivity_model == std::string("eucken")) &&
                     (viscosity_model == std::string("blottner")) )
              {
                physics_list[physics_to_add] = 
                  PhysicsPtr(new GRINS::ReactingLowMachNavierStokes<
                             GRINS::AntiochWilkeTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                 Antioch::MixtureViscosity<Antioch::BlottnerViscosity<libMesh::Real> >,
                                                                 Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                 Antioch::ConstantLewisDiffusivity<libMesh::Real> >,
                             GRINS::AntiochWilkeTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                   Antioch::MixtureViscosity<Antioch::BlottnerViscosity<libMesh::Real> >,
                                                                   Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                   Antioch::ConstantLewisDiffusivity<libMesh::Real> > >(physics_to_add,input) );
              }
            else
              {
                            std::cerr << "Error: Unknown Antioch model combination: "
                                      << "viscosity_model    = " << viscosity_model << std::endl
                                      << "conductivity_model = " << conductivity_model << std::endl
                                      << "diffusivity_model  = " << diffusivity_model << std::endl
                                      << "thermo_model       = " << thermo_model << std::endl;
                            libmesh_error();
              }
          }
        else if( mixing_model == std::string("constant") )
          {
            if( viscosity_model != std::string("constant") )
              {
                std::cerr << "Error: For constant mixing_model, viscosity model must be constant!"
                          << std::endl;
                libmesh_error();
              }

            if( diffusivity_model != std::string("constant_lewis") )
              {
                std::cerr << "Error: For constant mixing_model, diffusivity model must be constant_lewis!"
                          << std::endl;
                libmesh_error();
              }

            if( (thermo_model == std::string("stat_mech")) &&
                (conductivity_model == std::string("constant")) )
              {
                physics_list[physics_to_add] = 
                  PhysicsPtr(new GRINS::ReactingLowMachNavierStokes<GRINS::AntiochConstantTransportMixture<GRINS::ConstantConductivity>,
                                                                    GRINS::AntiochConstantTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>, GRINS::ConstantConductivity> >(physics_to_add,input) );
              }
            else if( (thermo_model == std::string("cea")) &&
                     (conductivity_model == std::string("constant")) )
              {
                physics_list[physics_to_add] = 
                  PhysicsPtr(new GRINS::ReactingLowMachNavierStokes<GRINS::AntiochConstantTransportMixture<GRINS::ConstantConductivity>,
                                                                    GRINS::AntiochConstantTransportEvaluator<Antioch::CEAEvaluator<libMesh::Real>, GRINS::ConstantConductivity> >(physics_to_add,input) );
              }
            else if( (thermo_model == std::string("stat_mech")) &&
                (conductivity_model == std::string("constant_prandtl")) )
              {
                physics_list[physics_to_add] = 
                  PhysicsPtr(new GRINS::ReactingLowMachNavierStokes<GRINS::AntiochConstantTransportMixture<GRINS::ConstantPrandtlConductivity>,
                                                                    GRINS::AntiochConstantTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>, GRINS::ConstantPrandtlConductivity> >(physics_to_add,input) );
              }
            else if( (thermo_model == std::string("cea")) &&
                     (conductivity_model == std::string("constant_prandtl")) )
              {
                physics_list[physics_to_add] = 
                  PhysicsPtr(new GRINS::ReactingLowMachNavierStokes<GRINS::AntiochConstantTransportMixture<GRINS::ConstantPrandtlConductivity>,
                                                                    GRINS::AntiochConstantTransportEvaluator<Antioch::CEAEvaluator<libMesh::Real>, GRINS::ConstantPrandtlConductivity> >(physics_to_add,input) );
              }
            else
              {
                std::cerr << "Error: Unknown Antioch model combination: "
                          << "viscosity_model    = " << viscosity_model << std::endl
                          << "conductivity_model = " << conductivity_model << std::endl
                          << "diffusivity_model  = " << diffusivity_model << std::endl
                          << "thermo_model       = " << thermo_model << std::endl;
                libmesh_error();
              }
          }
        else // mixing_model
          {
            std::cerr << "Error: Unknown Antioch mixing_model "
                      << mixing_model << "!" << std::endl;
            libmesh_error();
          }
#else
        std::cerr << "Error: Antioch not enabled. Cannot use Antioch library."
                  << std::endl;
        libmesh_error();

#endif // GRINS_HAVE_ANTIOCH
      }
    else
      {
        std::cerr << "Error: Invalid thermo-chemistry library" << std::endl
                  << "       for ReactingLowMachNavierStokes physics." << std::endl
                  << "       thermochemistry_library = " << thermochem_lib << std::endl;
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
        if( (physics->first == incompressible_navier_stokes_adjoint_stab) ||
            (physics->first == incompressible_navier_stokes_spgsm_stab) )
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

        /* For BoussinesqBuoyancyAdjointStabilization, we'd better have IncompressibleNavierStokes */
	if( physics->first == boussinesq_buoyancy_adjoint_stab )
	  {
	    if( physics_list.find(incompressible_navier_stokes) == physics_list.end() )
	      {
		this->physics_consistency_error( physics->first, incompressible_navier_stokes  );
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

  void PhysicsFactory::visc_error( const std::string& physics,				  
				   const std::string& viscosity ) const
  {
    std::cerr << "================================================================" << std::endl
	      << "Invalid combination of models for " << physics << std::endl	      
	      << "Viscosity model     = " << viscosity << std::endl	      
	      << "================================================================" << std::endl;
    libmesh_error();
  }

  void PhysicsFactory::conductivity_error( const std::string& physics,				  
				   const std::string& conductivity ) const
  {
    std::cerr << "================================================================" << std::endl
	      << "Invalid combination of models for " << physics << std::endl	      
	      << "Conductivity model     = " << conductivity << std::endl	      
	      << "================================================================" << std::endl;
    libmesh_error();
  }


} // namespace GRINS
