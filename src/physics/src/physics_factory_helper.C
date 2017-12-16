//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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

// These functions
#include "grins/physics_factory_helper.h"

// GRINS
#include "grins/common.h"
#include "grins/physics_naming.h"
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  void PhysicsFactoryHelper::parse_viscosity_model( const GetPot& input,
                                                    const std::string& physics,
                                                    std::string& model )
  {
    // Newer, preferred version
    bool have_material = MaterialsParsing::have_material(input,physics);

    // Old deprecated version
    bool have_viscosity_model = input.have_variable("Physics/"+PhysicsNaming::incompressible_navier_stokes()+"/viscosity_model");

    PhysicsFactoryHelper::deprecated_visc_model_parsing(  have_viscosity_model,
                                                          have_material,
                                                          input,
                                                          physics,
                                                          model );

    if( have_material )
      {
        std::string material = MaterialsParsing::material_name( input, physics );
        MaterialsParsing::viscosity_model( input, physics, material, model );
      }
  }

  void PhysicsFactoryHelper::parse_conductivity_model( const GetPot& input,
                                                       const std::string& physics,
                                                       std::string& model )
  {
    // Newer, preferred version
    bool have_material = MaterialsParsing::have_material(input,physics);

    // Old deprecated versions
    bool have_ht_conductivity_model = input.have_variable("Physics/"+PhysicsNaming::heat_transfer()+"/conductivity_model");

    bool have_conductivity_model = input.have_variable("Physics/"+physics+"/conductivity_model");

    if( (have_material && have_conductivity_model) ||
        (have_material && have_ht_conductivity_model) )
      {
        libmesh_error_msg("Error: Cannot specify both conductivity_model and material.");
      }

    /* If the user hasn't specified the material or the conductivity_model,
       we're assuming they're using the old version.
       This is deprecated.*/
    if( !have_conductivity_model && !have_material && !have_ht_conductivity_model )
      {
        std::string warning = "Warning: Neither conductivity_model nor material were specified.\n";
        warning += "      We are assuming a constant conductivity model.\n";
        warning += "      This case is DEPRECATED.\n";
        warning += "      Please update and specify Physics/"+physics+"/material.\n";
        grins_warning(warning);

        model = "constant";
      }

    // Deprecated
    if( have_ht_conductivity_model )
      {
        std::string warning = "Warning: Option Physics/"+PhysicsNaming::heat_transfer()+"/conductivity_model is DEPRECATED.\n";
        warning += "         Please update to use Physics/"+physics+"/material.\n";
        grins_warning(warning);

        model = input( "Physics/"+PhysicsNaming::heat_transfer()+"/conductivity_model", "constant" );
      }

    // Deprecated
    if( have_conductivity_model )
      {
        std::string warning = "Warning: Option Physics/"+physics+"/conductivity_model is DEPRECATED.\n";
        warning += "         Please update to use Physics/"+physics+"/material.\n";
        grins_warning(warning);

        model = input( "Physics/"+physics+"/conductivity_model", "constant" );
      }

    // Preferred
    if( have_material )
      {
        std::string material = MaterialsParsing::material_name( input, physics );
        MaterialsParsing::thermal_conductivity_model( input, physics, material, model );
      }

    return;
  }

  void PhysicsFactoryHelper::parse_specific_heat_model( const GetPot& input,
                                                        const std::string& physics,
                                                        std::string& model )
  {
    // Newer, preferred version
    bool have_material = MaterialsParsing::have_material(input,physics);

    bool have_lmns_specific_heat_model = input.have_variable("Physics/"+PhysicsNaming::low_mach_navier_stokes()+"/specific_heat_model");

    if( have_material && have_lmns_specific_heat_model )
      {
        libmesh_error_msg("ERROR: Cannot specify both a material and Physics/"+PhysicsNaming::low_mach_navier_stokes()+"/specific_heat_model!");
      }

    // Deprecated
    if( have_lmns_specific_heat_model )
      {
        std::string warning = "Warning: Option Physics/"+PhysicsNaming::low_mach_navier_stokes()+"/specific_heat_model is DEPRECATED.\n";
        warning += "         Please update to use Physics/"+physics+"/material.\n";
        grins_warning(warning);

        model = input( "Physics/"+PhysicsNaming::low_mach_navier_stokes()+"/specific_heat_model", "constant" );
      }

    // Preferred
    if( have_material )
      {
        std::string material = MaterialsParsing::material_name( input, physics );
        MaterialsParsing::specific_heat_model( input, physics, material, model );
      }
  }

  void PhysicsFactoryHelper::parse_turb_viscosity_model( const GetPot& input,
                                                         const std::string& physics,
                                                         std::string& model )
  {
    // Newer, preferred version
    bool have_material = MaterialsParsing::have_material(input,physics);

    // Old deprecated version
    bool have_viscosity_model = input.have_variable("Physics/"+PhysicsNaming::incompressible_navier_stokes()+"/viscosity_model");

    PhysicsFactoryHelper::deprecated_visc_model_parsing(  have_viscosity_model,
                                                          have_material,
                                                          input,
                                                          physics,
                                                          model );

    // Additionally, we need to check if the turbulence Physics didn't
    // set the viscosity model, but the *non* turbulent viscosity_model
    // is set to spallartallmaras. In that case, the physical viscosity
    // model is 'constant'.
    if( !have_material &&
        input.have_variable( "Physics/"+PhysicsNaming::incompressible_navier_stokes()+"/viscosity_model") )
      {
        // If we got here, have_viscosity_model is true and might
        // be set to spallartallmaras
        if( input("Physics/"+PhysicsNaming::incompressible_navier_stokes()+"/viscosity_model", "DIE!") == std::string("spalartallmaras") )
          {
            model = "constant";
          }
      }

    if( have_material )
      {
        std::string material = MaterialsParsing::material_name( input, physics );
        MaterialsParsing::turb_viscosity_model( input, physics, material, model );
      }

    return;
  }

  void PhysicsFactoryHelper::parse_stress_strain_model( const GetPot& input,
                                                        const std::string& physics,
                                                        std::string& model,
                                                        std::string& strain_energy )
  {
    // Newer, preferred version
    std::string material = MaterialsParsing::material_name( input, physics );

    // Old deprecated version
    bool have_elasticity_model = input.have_variable("Physics/"+physics+"/elasticity_model");

    // It's an error to specify both the old and the new version
    if( have_elasticity_model &&
        input.have_variable("Materials/"+material+"/StressStrainLaw/model") )
      {
        libmesh_error_msg("ERROR: Cannot specify both Materials/"+material+"/StressStrainLaw/model and Physics/"+physics+"/elasticity_model!");
      }

    // It's an error if don't specify at least one of them
    if( !have_elasticity_model &&
        !input.have_variable("Materials/"+material+"/StressStrainLaw/model") )
      {
        // But since the old is deprecated, we'll just them to supply the new
        libmesh_error_msg("ERROR: Must specify Materials/"+material+"/StressStrainLaw/model!");
      }

    // Deprecated
    if( have_elasticity_model )
      {
        std::string warning = "Warning: Option Physics/"+physics+"/elasticity_model is DEPRECATED.\n";
        warning += "         Please update to use Materials/MATERIAL_NAME/StressStrainLaw/model.\n";
        grins_warning(warning);

        model = input( "Physics/"+physics+"/elasticity_model", "DIE!" );

        if( model == std::string("HookesLaw") )
          model = "hookes_law";
        if( model == std::string("MooneyRivlin") )
          {
            model = "incompressible_hyperelasticity";
            strain_energy = "mooney_rivlin";
          }

      }
    // Preferred
    else if( input.have_variable("Materials/"+material+"/StressStrainLaw/model") )
      {
        MaterialsParsing::stress_strain_model( input, physics, material,
                                               model, strain_energy );
      }
    // Wat?
    else
      {
        libmesh_error();
      }
  }

  void PhysicsFactoryHelper::parse_thermochemistry_model( const GetPot& input,
                                                          const std::string& physics,
                                                          std::string& model )
  {
    // Newer, preferred version
    std::string material = MaterialsParsing::material_name( input, physics );

    bool have_thermochem_lib = input.have_variable( "Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/thermochemistry_library" );

    // It's an error to specify both the old and the new version
    if( have_thermochem_lib &&
        input.have_variable("Materials/"+material+"/GasMixture/thermochemistry_library") )
      {
        libmesh_error_msg("ERROR: Cannot specify both Materials/"+material+"/GasMixture/thermochemistry_library and Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/thermochemistry_library!");
      }

    //Deprecated
    if( have_thermochem_lib )
      {
        model = input( "Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/thermochemistry_library", "DIE!" );
      }
    // Preferred
    else if( input.have_variable("Materials/"+material+"/GasMixture/thermochemistry_library") )
      {
        model = input("Materials/"+material+"/GasMixture/thermochemistry_library", "DIE!");
      }
    // Fail
    else
      {
        libmesh_error_msg("ERROR! Could not find valid thermochemistry_library input!");
      }

    // Make sure we have a valid model
    if( model != std::string("antioch") &&
        model != std::string("cantera") )
      {
        std::string error = "ERROR! Invalid thermochemistry library value "+model+"!\n";
        error += "       Valid selections are: antioch\n";
        error += "                             cantera\n";

        libmesh_error_msg(error);
      }
  }

  void PhysicsFactoryHelper::parse_antioch_models( const GetPot& input,
                                                   const std::string& physics,
                                                   std::string& transport_model,
                                                   std::string& thermo_model,
                                                   std::string& viscosity_model,
                                                   std::string& conductivity_model,
                                                   std::string& diffusivity_model )
  {
    // Newer, preferred version
    std::string material = MaterialsParsing::material_name( input, physics );

    bool have_transport_model = input.have_variable( "Physics/Antioch/transport_model" );

    // It's an error to specify both the old and the new version
    if( have_transport_model &&
        input.have_variable("Materials/"+material+"/GasMixture/Antioch/transport_model") )
      {
        libmesh_error_msg("ERROR: Cannot specify both Materials/"+material+"/GasMixture/Antioch/transport_model and Physics/Antioch/transport_model!");
      }

    //Deprecated
    if( have_transport_model )
      {
        std::string warning = "Warning: Option Physics/Antioch/transport_model is DEPRECATED.\n";
        warning += "         Please update to use Use Materials/MATERIAL_NAME/GasMixture/Antioch/transport_model.\n";
        grins_warning(warning);

        transport_model = input( "Physics/Antioch/transport_model", "mixture_averaged" );
      }
    // mixing_model option is now deprecated in favor of transport_model
    else if( input.have_variable("Physics/Antioch/mixing_model") )
      {
        std::string warning = "Warning: Option Physics/Antioch/mixing_model is DEPRECATED.\n";
        warning += "         Please update to use Use Materials/MATERIAL_NAME/GasMixture/Antioch/transport_model.\n";
        grins_warning(warning);

        transport_model = input( "Physics/Antioch/mixing_model" , "mixture_averaged" );
      }
    // Preferred
    else if( input.have_variable("Materials/"+material+"/GasMixture/Antioch/transport_model") )
      {
        transport_model = input("Materials/"+material+"/GasMixture/Antioch/transport_model", "DIE!");
      }
    // Fail
    else
      {
        libmesh_error_msg("ERROR! Could not find valid transport_model input!");
      }

    // transport_model = wilke is deprecated
    if( transport_model == std::string("wilke") )
      {
        libMesh::err << "WARNING: Physics/Antioch/transport_model value of 'wilke' is deprecated!" << std::endl
                     << "         Replace Physics/Antioch/transport_model value with 'mixture_averaged'"
                     << std::endl;

        transport_model = "mixture_averaged";
      }

    // Now parse the remaining models
    //Deprecated
    if( have_transport_model || input.have_variable("Physics/Antioch/mixing_model") )
      {
        // We're tying all the other option specifications to transport_model.
        thermo_model = input( "Physics/Antioch/thermo_model", "stat_mech");
        viscosity_model = input( "Physics/Antioch/viscosity_model", "blottner");
        conductivity_model = input( "Physics/Antioch/conductivity_model", "eucken");
        diffusivity_model = input( "Physics/Antioch/diffusivity_model", "constant_lewis");

        // So they'd better not have the other version specified.
        if( input.have_variable("Materials/"+material+"/GasMixture/Antioch/thermo_model")       ||
            input.have_variable("Materials/"+material+"/GasMixture/Antioch/viscosity_model")    ||
            input.have_variable("Materials/"+material+"/GasMixture/Antioch/thermal_conductivity_model") ||
            input.have_variable("Materials/"+material+"/GasMixture/Antioch/mass_diffusivity_model") )
          {
            libmesh_error_msg("ERROR: Cannot specifiy Physics/Antioch/transport_model and then specify Materials/"+material+"/GasMixture/Antioch/<thermo,viscosity,conductivity,diffusivity>_model!");
          }
      }
    // Preferred
    else if( input.have_variable("Materials/"+material+"/GasMixture/Antioch/transport_model") )
      {
        // We're tying all the other option specifications to transport_model.
        // For the newer cases, we don't support a default, the user must specify
        // We're tying all the other option specifications to transport_model.
        MaterialsParsing::check_for_input_option(input,"Materials/"+material+"/GasMixture/Antioch/thermo_model");
        MaterialsParsing::check_for_input_option(input,"Materials/"+material+"/GasMixture/Antioch/viscosity_model");
        MaterialsParsing::check_for_input_option(input,"Materials/"+material+"/GasMixture/Antioch/thermal_conductivity_model");
        MaterialsParsing::check_for_input_option(input,"Materials/"+material+"/GasMixture/Antioch/mass_diffusivity_model");

        thermo_model = input( "Materials/"+material+"/GasMixture/Antioch/thermo_model", "DIE!");
        viscosity_model = input( "Materials/"+material+"/GasMixture/Antioch/viscosity_model", "DIE!");
        conductivity_model = input( "Materials/"+material+"/GasMixture/Antioch/thermal_conductivity_model", "DIE!");
        diffusivity_model = input( "Materials/"+material+"/GasMixture/Antioch/mass_diffusivity_model", "DIE!");

        // So they'd better not have the other version specified.
        if( input.have_variable("Physics/Antioch/thermo_model")       ||
            input.have_variable("Physics/Antioch/viscosity_model")    ||
            input.have_variable("Physics/Antioch/conductivity_model") ||
            input.have_variable("Physics/Antioch/diffusivity_model") )
          {
            libmesh_error_msg("ERROR: Cannot specifiy Materials/"+material+"/GasMixture/Antioch/transport_model and then specify Physics/Antioch/<thermo,viscosity,conductivity,diffusivity>_model!");
          }
      }
    // Fail
    else
      {
        libmesh_error_msg("ERROR! Could not find valid transport_model input!");
      }

  }

  void PhysicsFactoryHelper::deprecated_visc_model_parsing( bool have_ins_viscosity_model,
                                                            bool have_material,
                                                            const GetPot& input,
                                                            const std::string& physics,
                                                            std::string& model)
  {
    bool have_viscosity_model = input.have_variable( "Physics/"+physics+"/viscosity_model");

    if( (have_material && have_ins_viscosity_model) ||
        (have_material && have_viscosity_model) )
      {
        libmesh_error_msg("Error: Cannot specify both viscosity_model and material.");
      }

    /* If the user hasn't specified the material or the viscosity_model,
       we're assuming they're using the old version. */
    if( !have_ins_viscosity_model && !have_material && !have_viscosity_model )
      {
        std::string warning = "Warning: Neither viscosity_model nor material were specified.\n";
        warning += "      We are assuming a constant viscosity model.\n";
        warning += "      This case is DEPRECATED.\n";
        warning += "      Please update and specify Physics/"+physics+"/material.\n";
        grins_warning(warning);

        model = "constant";
      }

    if( have_ins_viscosity_model )
      {
        std::string warning = "Warning: Option Physics/"+PhysicsNaming::incompressible_navier_stokes()+"/viscosity_model is DEPRECATED.\n";
        warning += "         Please update to use Physics/"+physics+"/material.\n";
        grins_warning(warning);

        model = input( "Physics/"+PhysicsNaming::incompressible_navier_stokes()+"/viscosity_model", "constant" );
      }

    if( have_viscosity_model )
      {
        std::string warning = "Warning: Option Physics/"+physics+"/viscosity_model is DEPRECATED.\n";
        warning += "         Please update to use Physics/"+physics+"/material.\n";
        grins_warning(warning);

        model = input( "Physics/"+physics+"/viscosity_model", "constant" );
      }
  }

} // end namespace GRINS
