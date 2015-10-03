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


// This class
#include "grins/physics_factory.h"
#include "grins/physics_factory_helper.h"

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
#include "grins/averaged_fan_adjoint_stab.h"
#include "grins/averaged_turbine.h"
#include "grins/scalar_ode.h"
#include "grins/velocity_drag.h"
#include "grins/velocity_drag_adjoint_stab.h"
#include "grins/velocity_penalty.h"
#include "grins/velocity_penalty_adjoint_stab.h"
#include "grins/parsed_velocity_source.h"
#include "grins/parsed_velocity_source_adjoint_stab.h"
#include "grins/elastic_membrane.h"
#include "grins/elastic_cable.h"
#include "grins/elastic_membrane_constant_pressure.h"
#include "grins/elastic_cable_constant_gravity.h"
#include "grins/grins_physics_names.h"
#include "grins/constant_source_term.h"
#include "grins/parsed_source_term.h"

#include "grins/spalart_allmaras.h"
#include "grins/spalart_allmaras_spgsm_stab.h"

#include "grins/constant_conductivity.h"
#include "grins/constant_specific_heat.h"
#include "grins/constant_viscosity.h"
#include "grins/parsed_conductivity.h"
#include "grins/parsed_viscosity.h"
#include "grins/spalart_allmaras_viscosity.h"

#include "grins/reacting_low_mach_navier_stokes.h"
#include "grins/heat_conduction.h"
#include "grins/constant_source_func.h"

#include "grins/antioch_mixture_averaged_transport_evaluator.h"
#include "grins/antioch_constant_transport_mixture.h"
#include "grins/antioch_constant_transport_evaluator.h"

#include "grins/hookes_law.h"
#include "grins/hookes_law_1d.h"
#include "grins/incompressible_plane_stress_hyperelasticity.h"
#include "grins/mooney_rivlin.h"

// Antioch
#ifdef GRINS_HAVE_ANTIOCH
#include "antioch_config.h"
#include "antioch/sutherland_viscosity.h"
#include "antioch/blottner_viscosity.h"
#include "antioch/sutherland_parsing.h"
#include "antioch/blottner_parsing.h"
#include "antioch/eucken_thermal_conductivity.h"
#include "antioch/constant_lewis_diffusivity.h"
#endif // GRINS_HAVE_ANTIOCH

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{

  // Starting with helper functions for model-templated physics
  // instantiation.  These don't need to be visible elsewhere.

  typedef PhysicsFactory::PhysicsPtr PhysicsPtr;

  void visc_cond_specheat_error( const std::string& physics,
				 const std::string& conductivity,
				 const std::string& viscosity,
				 const std::string& specific_heat )
  {
    std::cerr << "================================================================" << std::endl
	      << "Invalid combination of models for " << physics << std::endl
	      << "Conductivity model  = " << conductivity << std::endl
	      << "Viscosity model     = " << viscosity << std::endl
	      << "Specific heat model = " << specific_heat << std::endl
	      << "================================================================" << std::endl;
    libmesh_error();
  }

  void visc_error( const std::string& physics,
		   const std::string& viscosity )
  {
    std::cerr << "================================================================" << std::endl
	      << "Invalid combination of models for " << physics << std::endl
	      << "Viscosity model     = " << viscosity << std::endl
	      << "================================================================" << std::endl;
    libmesh_error();
  }

  void conductivity_error( const std::string& physics,
			   const std::string& conductivity )
  {
    std::cerr << "================================================================" << std::endl
	      << "Invalid combination of models for " << physics << std::endl
	      << "Conductivity model     = " << conductivity << std::endl
	      << "================================================================" << std::endl;
    libmesh_error();
  }


  // FIXME: Currently we always look for a lone viscosity model in the
  // IncompressibleNavierStokes input, or for a lone conductivity
  // model in the HeatTransfer input, or for a viscosity + specific
  // heat + conductivity model set in the LowMachNavierStokes input.
  // This needs to be fixed, but carefully to avoid breaking old input
  // files.

  // We use the core_physics to look for the material model
  template <template<typename> class Subclass>
  PhysicsPtr new_mu_class(const std::string& physics_to_add,
                          const std::string& core_physics,
                          const GetPot& input)
  {
    std::string viscosity;
    PhysicsFactoryHelper::parse_viscosity_model(input,core_physics,viscosity);

    if( viscosity == "constant" )
      return PhysicsPtr
        (new Subclass<ConstantViscosity>(physics_to_add,input));
    else if( viscosity == "parsed" )
      return PhysicsPtr
        (new Subclass<ParsedViscosity>(physics_to_add,input));
    // For SA viscosity model, we need to parse what the "sub" viscosity model is
    else if( viscosity == "spalartallmaras" )
      {
        std::string turb_viscosity;
        PhysicsFactoryHelper::parse_turb_viscosity_model(input,core_physics,turb_viscosity);
        if( turb_viscosity == "constant" )
          return PhysicsPtr
            (new Subclass<SpalartAllmarasViscosity<ConstantViscosity> >(physics_to_add,input));
        visc_error(physics_to_add, turb_viscosity);
      }

    visc_error(physics_to_add, viscosity);
    return PhysicsPtr();
  }

  // Parse viscosity model for turbulence classes
  template <template<typename> class Subclass>
  PhysicsPtr new_turb_mu_class(const std::string& physics_to_add,
                               const std::string& core_physics,
                               const GetPot& input)
  {
    std::string viscosity;
    PhysicsFactoryHelper::parse_turb_viscosity_model(input,core_physics,viscosity);

    if( viscosity == "constant" )
      return PhysicsPtr
        (new Subclass<ConstantViscosity>(physics_to_add,input));

    visc_error(physics_to_add, viscosity);
    return PhysicsPtr();
  }

  // core_physics is used to determine material model
  template <template<typename> class Subclass>
  PhysicsPtr new_k_class(const std::string& physics_to_add,
                         const std::string& core_physics,
                         const GetPot& input)
  {
    std::string conductivity;
    PhysicsFactoryHelper::parse_conductivity_model(input,core_physics,conductivity);

    if( conductivity == "constant" )
      return PhysicsPtr
        (new Subclass<ConstantConductivity>(physics_to_add,input));
    else if( conductivity == "parsed" )
      return PhysicsPtr
        (new Subclass<ParsedConductivity>(physics_to_add,input));

    conductivity_error(physics_to_add, conductivity);
    return PhysicsPtr();
  }

  template <template<typename,typename,typename> class Subclass>
  PhysicsPtr new_mu_cp_k_class(const std::string& physics_to_add,
                               const std::string& core_physics,
                               const GetPot& input)
  {
    std::string conductivity;
    PhysicsFactoryHelper::parse_conductivity_model(input,core_physics,conductivity);

    std::string viscosity;
    PhysicsFactoryHelper::parse_viscosity_model(input,core_physics,viscosity);

    std::string specific_heat;
    PhysicsFactoryHelper::parse_specific_heat_model(input,core_physics,specific_heat);

    if(  conductivity == "constant" && viscosity == "constant" && specific_heat == "constant" )
      {
        return PhysicsPtr
          (new Subclass<ConstantViscosity,
                        ConstantSpecificHeat,
                        ConstantConductivity>
             (physics_to_add,input));
      }

    visc_cond_specheat_error(physics_to_add, conductivity,
                             viscosity, specific_heat);
    return PhysicsPtr();
  }

  template <template<typename> class Subclass>
  PhysicsPtr new_plane_stress_class(const std::string& physics_to_add,
                                    const std::string& core_physics,
                                    const GetPot& input)
  {
    std::string model = "none";
    std::string strain_energy = "none";

    PhysicsFactoryHelper::parse_stress_strain_model( input,
                                                     core_physics,
                                                     model,
                                                     strain_energy );

    if( model == std::string("hookes_law") )
      {
        return PhysicsPtr
          (new Subclass<HookesLaw>
           (physics_to_add,input, false /*is_compressible*/));
      }
    else if( model == std::string("incompressible_hyperelasticity") )
      {
        if( strain_energy == std::string("mooney_rivlin") )
          {
            return PhysicsPtr
              (new Subclass<IncompressiblePlaneStressHyperelasticity<MooneyRivlin> >
               (physics_to_add,input,false /*is_compressible*/));
          }
        else
          {
            std::string error = "ERROR: Invalid strain_energy "+strain_energy+"!\n";
            error += "       Valid values are: mooney_rivlin\n";
            libmesh_error_msg(error);
          }

      }
    else
      {
        std::string error = "Error: Invalid stress-strain model: "+model+"!\n";
        error += "       Valid values are: hookes_law\n";
        error += "                         incompressible_hyperelasticity\n";
        libmesh_error_msg(error);
      }

    // dummy
    return PhysicsPtr();
  }

  template <template<typename,typename> class Subclass>
  PhysicsPtr new_reacting_low_mach_class(const std::string& physics_to_add,
                                         const GetPot& input)
  {
    std::string thermochem_lib = input( "Physics/"+reacting_low_mach_navier_stokes+"/thermochemistry_library", "DIE!" );

    if( thermochem_lib == "cantera" )
      {
#ifdef GRINS_HAVE_CANTERA
        return PhysicsPtr(new Subclass<CanteraMixture,CanteraEvaluator>(physics_to_add,input));
#else
        std::cerr << "Error: Cantera not enabled. Cannot use Cantera library."
                  << std::endl;
        libmesh_error();

#endif // GRINS_HAVE_CANTERA
      }
    else if( thermochem_lib == "antioch" )
      {
#ifdef GRINS_HAVE_ANTIOCH

        std::string transport_model = input( "Physics/Antioch/transport_model" , "mixture_averaged" );

        // mixing_model option is now deprecated in favor of transport_model
        if( input.have_variable("Physics/Antioch/mixing_model") )
          {
            libMesh::err << "WARNING: Option Physics/Antioch/mixing_model is deprecated!" << std::endl
                         << "         Use Physics/Antioch/transport_model instead!" << std::endl;

            transport_model = input( "Physics/Antioch/mixing_model" , "mixture_averaged" );
          }

        // transport_model = wilke is deprecated
        if( transport_model == std::string("wilke") )
          {
            libMesh::err << "WARNING: Physics/Antioch/transport_model value of 'wilke' is deprecated!" << std::endl
                         << "         Replace Physics/Antioch/transport_model value with 'mixture_averaged'"
                         << std::endl;

            transport_model = "mixture_averaged";
          }

        std::string thermo_model = input( "Physics/Antioch/thermo_model", "stat_mech");
        std::string viscosity_model = input( "Physics/Antioch/viscosity_model", "blottner");
        std::string conductivity_model = input( "Physics/Antioch/conductivity_model", "eucken");
        std::string diffusivity_model = input( "Physics/Antioch/diffusivity_model", "constant_lewis");

        if( transport_model == std::string("mixture_averaged") )
          {
            if( (thermo_model == std::string("stat_mech")) &&
                (diffusivity_model == std::string("constant_lewis")) &&
                (conductivity_model == std::string("eucken")) &&
                (viscosity_model == std::string("sutherland")) )
              {
                return PhysicsPtr(new Subclass<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                             Antioch::SutherlandViscosity<libMesh::Real>,
                                                                                             Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                                             Antioch::ConstantLewisDiffusivity<libMesh::Real> >,
                                               GRINS::AntiochMixtureAveragedTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                               Antioch::SutherlandViscosity<libMesh::Real>,
                                                                                               Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                                               Antioch::ConstantLewisDiffusivity<libMesh::Real> > >(physics_to_add,input) );
              }
            else if( (thermo_model == std::string("stat_mech")) &&
                     (diffusivity_model == std::string("constant_lewis")) &&
                     (conductivity_model == std::string("eucken")) &&
                     (viscosity_model == std::string("blottner")) )
              {
                return PhysicsPtr(new Subclass<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                             Antioch::BlottnerViscosity<libMesh::Real>,
                                                                                             Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                                             Antioch::ConstantLewisDiffusivity<libMesh::Real> >,
                                               GRINS::AntiochMixtureAveragedTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                               Antioch::BlottnerViscosity<libMesh::Real>,
                                                                                               Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                                               Antioch::ConstantLewisDiffusivity<libMesh::Real> > >(physics_to_add,input) );
              }
            else if( (thermo_model == std::string("stat_mech")) &&
                     (diffusivity_model == std::string("kinetics_theory")) &&
                     (conductivity_model == std::string("kinetics_theory")) &&
                     (viscosity_model == std::string("kinetics_theory")) )
              {
#ifdef ANTIOCH_HAVE_GSL
                return PhysicsPtr(new Subclass<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                   Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>,
                                                                                   Antioch::KineticsTheoryThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real>,libMesh::Real>,
                                                                                   Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> >,
                                               GRINS::AntiochMixtureAveragedTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                     Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>,
                                                                                     Antioch::KineticsTheoryThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real>,libMesh::Real>,
                                                                                     Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> > >(physics_to_add,input) );
#else
                libmesh_error_msg("ERROR: Antioch requires GSL in order to use kinetics theory based models!");
#endif // ANTIOCH_HAVE_GSL
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
        else if( transport_model == std::string("constant") )
          {
            if( viscosity_model != std::string("constant") )
              {
                std::cerr << "Error: For constant transport_model, viscosity model must be constant!"
                          << std::endl;
                libmesh_error();
              }

            if( diffusivity_model != std::string("constant_lewis") )
              {
                std::cerr << "Error: For constant transport_model, diffusivity model must be constant_lewis!"
                          << std::endl;
                libmesh_error();
              }

            if( (thermo_model == std::string("stat_mech")) &&
                (conductivity_model == std::string("constant")) )
              {
                return PhysicsPtr(new Subclass<GRINS::AntiochConstantTransportMixture<GRINS::ConstantConductivity>,
                                               GRINS::AntiochConstantTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>, GRINS::ConstantConductivity> >(physics_to_add,input) );
              }
            else if( (thermo_model == std::string("cea")) &&
                     (conductivity_model == std::string("constant")) )
              {
                return PhysicsPtr(new Subclass<GRINS::AntiochConstantTransportMixture<GRINS::ConstantConductivity>,
                                               GRINS::AntiochConstantTransportEvaluator<Antioch::CEAEvaluator<libMesh::Real>, GRINS::ConstantConductivity> >(physics_to_add,input) );
              }
            else if( (thermo_model == std::string("stat_mech")) &&
                (conductivity_model == std::string("constant_prandtl")) )
              {
                return PhysicsPtr(new Subclass<GRINS::AntiochConstantTransportMixture<GRINS::ConstantPrandtlConductivity>,
                                               GRINS::AntiochConstantTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>, GRINS::ConstantPrandtlConductivity> >(physics_to_add,input) );
              }
            else if( (thermo_model == std::string("cea")) &&
                     (conductivity_model == std::string("constant_prandtl")) )
              {
                return PhysicsPtr(new Subclass<GRINS::AntiochConstantTransportMixture<GRINS::ConstantPrandtlConductivity>,
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
        else // transport_model
          {
            std::cerr << "Error: Unknown Antioch transport_model "
                      << transport_model << "!" << std::endl
                      << "       Valid values are: constant" << std::endl
                      << "                         mixture_averaged" << std::endl;
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

    return PhysicsPtr();
  }


  // And now PhysicsFactory methods

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
	physics_list[physics_to_add] =
          new_mu_class<IncompressibleNavierStokes>
          (physics_to_add, incompressible_navier_stokes, input);
      }
    else if( physics_to_add == stokes )
      {
	physics_list[physics_to_add] =
          new_mu_class<Stokes>
          (physics_to_add, stokes, input);
      }
    else if( physics_to_add == incompressible_navier_stokes_adjoint_stab )
      {
	physics_list[physics_to_add] =
          new_mu_class<IncompressibleNavierStokesAdjointStabilization>
          (physics_to_add, incompressible_navier_stokes, input);
      }
    else if( physics_to_add == incompressible_navier_stokes_spgsm_stab )
      {
	physics_list[physics_to_add] =
          new_mu_class<IncompressibleNavierStokesSPGSMStabilization>
          (physics_to_add, incompressible_navier_stokes, input);
      }
    else if( physics_to_add == velocity_drag )
      {
	physics_list[physics_to_add] =
          new_mu_class<VelocityDrag>
          (physics_to_add, incompressible_navier_stokes, input);
      }
    else if( physics_to_add == velocity_drag_adjoint_stab )
      {
	physics_list[physics_to_add] =
          new_mu_class<VelocityDragAdjointStabilization>
          (physics_to_add, incompressible_navier_stokes, input);
      }
    else if( physics_to_add == velocity_penalty  ||
             physics_to_add == velocity_penalty2 ||
             physics_to_add == velocity_penalty3)
      {
	physics_list[physics_to_add] =
          new_mu_class<VelocityPenalty>
          (physics_to_add, incompressible_navier_stokes, input);
      }
    else if( physics_to_add == velocity_penalty_adjoint_stab  ||
             physics_to_add == velocity_penalty2_adjoint_stab ||
             physics_to_add == velocity_penalty3_adjoint_stab )
      {
	physics_list[physics_to_add] =
          new_mu_class<VelocityPenaltyAdjointStabilization>
          (physics_to_add, incompressible_navier_stokes, input);
      }
    else if( physics_to_add == parsed_velocity_source )
      {
	physics_list[physics_to_add] =
          new_mu_class<ParsedVelocitySource>
          (physics_to_add, incompressible_navier_stokes, input);
      }
    else if( physics_to_add == parsed_velocity_source_adjoint_stab )
      {
	physics_list[physics_to_add] =
          new_mu_class<ParsedVelocitySourceAdjointStabilization>
          (physics_to_add, incompressible_navier_stokes, input);
      }
    else if( physics_to_add == averaged_fan )
      {
	physics_list[physics_to_add] =
          new_mu_class<AveragedFan>
          (physics_to_add, incompressible_navier_stokes, input);
      }
    else if( physics_to_add == averaged_fan_adjoint_stab )
      {
	physics_list[physics_to_add] =
          new_mu_class<AveragedFanAdjointStabilization>
          (physics_to_add, incompressible_navier_stokes, input);
      }
    else if( physics_to_add == averaged_turbine )
      {
	physics_list[physics_to_add] =
          new_mu_class<AveragedTurbine>
          (physics_to_add, incompressible_navier_stokes, input);
      }
    else if( physics_to_add == spalart_allmaras )
      {
        physics_list[physics_to_add] =
          new_turb_mu_class<SpalartAllmaras>
          (physics_to_add, spalart_allmaras, input);
      }
    else if( physics_to_add == spalart_allmaras_spgsm_stab )
      {
        physics_list[physics_to_add] =
          new_turb_mu_class<SpalartAllmarasSPGSMStabilization>
          (physics_to_add, spalart_allmaras, input);
      }
    else if( physics_to_add == scalar_ode )
      {
	physics_list[physics_to_add] =
	  PhysicsPtr(new ScalarODE(physics_to_add,input));
      }
    else if( physics_to_add == heat_transfer )
      {
	physics_list[physics_to_add] =
          new_k_class<HeatTransfer>
          (physics_to_add, heat_transfer, input);
      }
    else if( physics_to_add == heat_transfer_adjoint_stab )
      {
	physics_list[physics_to_add] =
          new_k_class<HeatTransferAdjointStabilization>
          (physics_to_add, heat_transfer, input);
      }
    else if( physics_to_add == heat_transfer_spgsm_stab )
      {
	physics_list[physics_to_add] =
          new_k_class<HeatTransferSPGSMStabilization>
          (physics_to_add, heat_transfer, input);
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
	physics_list[physics_to_add] =
          new_k_class<AxisymmetricHeatTransfer>
          (physics_to_add, axisymmetric_heat_transfer, input);
      }
    else if( physics_to_add == boussinesq_buoyancy )
      {
	physics_list[physics_to_add] =
	  PhysicsPtr(new BoussinesqBuoyancy(physics_to_add,input));
      }
    else if( physics_to_add == boussinesq_buoyancy_adjoint_stab )
      {
        physics_list[physics_to_add] =
          new_mu_class<BoussinesqBuoyancyAdjointStabilization>
          (physics_to_add, boussinesq_buoyancy, input);
      }
    else if( physics_to_add == boussinesq_buoyancy_spgsm_stab )
      {
        physics_list[physics_to_add] =
          new_mu_class<BoussinesqBuoyancySPGSMStabilization>
          (physics_to_add, boussinesq_buoyancy, input);
      }
    else if( physics_to_add == axisymmetric_boussinesq_buoyancy)
      {
	physics_list[physics_to_add] =
	  PhysicsPtr(new AxisymmetricBoussinesqBuoyancy(physics_to_add,input));
      }
    else if( physics_to_add == heat_conduction )
      {
	physics_list[physics_to_add] =
          new_k_class<HeatConduction>
          (physics_to_add, heat_conduction, input);
      }
    else if(  physics_to_add == low_mach_navier_stokes )
      {
	physics_list[physics_to_add] =
          new_mu_cp_k_class<LowMachNavierStokes>
          (physics_to_add, low_mach_navier_stokes, input);
      }
    else if(  physics_to_add == low_mach_navier_stokes_spgsm_stab )
      {
	physics_list[physics_to_add] =
          new_mu_cp_k_class<LowMachNavierStokesSPGSMStabilization>
          (physics_to_add, low_mach_navier_stokes, input);
      }
    else if(  physics_to_add == low_mach_navier_stokes_vms_stab )
      {
	physics_list[physics_to_add] =
          new_mu_cp_k_class<LowMachNavierStokesVMSStabilization>
            (physics_to_add, low_mach_navier_stokes, input);
      }
    else if(  physics_to_add == low_mach_navier_stokes_braack_stab )
      {
	physics_list[physics_to_add] =
          new_mu_cp_k_class<LowMachNavierStokesBraackStabilization>
            (physics_to_add, low_mach_navier_stokes, input);
      }
    else if( physics_to_add == reacting_low_mach_navier_stokes )
      {
        physics_list[physics_to_add] =
          new_reacting_low_mach_class<ReactingLowMachNavierStokes>
             (physics_to_add, input);
      }
    else if( physics_to_add == elastic_membrane )
      {
        physics_list[physics_to_add] =
          new_plane_stress_class<ElasticMembrane>
          (physics_to_add, elastic_membrane, input);
      }
    else if( physics_to_add == elastic_membrane_constant_pressure )
      {
        physics_list[physics_to_add] =
          PhysicsPtr(new ElasticMembraneConstantPressure(physics_to_add,input));
      }
    else if( physics_to_add == elastic_cable )
      {
        std::string elasticity_model = input("Physics/"+elastic_cable+"/elasticity_model", "HookesLaw" );

        if( elasticity_model == std::string("HookesLaw") )
          {
            physics_list[physics_to_add] =
              PhysicsPtr(new ElasticCable<HookesLaw1D>(physics_to_add,input,false /*is_compressible*/));
          }
        else
          {
            std::cerr << "Error: Invalid elasticity_model: " << elasticity_model << std::endl
                      << "       Valid selections are: Hookean" << std::endl;
            libmesh_error();
          }
      }
    else if( physics_to_add == elastic_cable_constant_gravity )
      {
        physics_list[physics_to_add] =
          PhysicsPtr(new ElasticCableConstantGravity(physics_to_add,input));
      }
    else if( physics_to_add == constant_source_term )
      {
        physics_list[physics_to_add] =
          PhysicsPtr(new ConstantSourceTerm(physics_to_add,input));
      }
    else if( physics_to_add == parsed_source_term )
      {
        physics_list[physics_to_add] =
          PhysicsPtr(new ParsedSourceTerm(physics_to_add,input));
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

} // namespace GRINS
