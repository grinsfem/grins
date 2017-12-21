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
#include "grins/materials_parsing.h"

// GRINS
#include "grins/common.h"
#include "grins/physics_naming.h"
#include "grins/parameter_user.h"

namespace GRINS
{
  void MaterialsParsing::viscosity_model( const GetPot & input,
                                          const std::string & physics,
                                          std::string & model )
  {
    std::string material = MaterialsParsing::material_name( input, physics );
    std::string visc_option("Materials/"+material+"/Viscosity/model");

    MaterialsParsing::check_for_input_option(input,visc_option);

    model = input( visc_option, "DIE!" );
  }

  void MaterialsParsing::thermal_conductivity_model( const GetPot & input,
                                                     const std::string & physics,
                                                     std::string & model )
  {
    std::string material = MaterialsParsing::material_name( input, physics );
    std::string cond_option("Materials/"+material+"/ThermalConductivity/model");

    MaterialsParsing::check_for_input_option(input,cond_option);

    model = input( cond_option, "DIE!" );
  }

  void MaterialsParsing::specific_heat_model( const GetPot & input,
                                              const std::string & /*physics*/,
                                              const std::string & material,
                                              std::string & model )
  {
    if( !input.have_variable("Materials/"+material+"/SpecificHeat/model") )
      {
        libmesh_error_msg("Error: Could not find Materials/"+material+"/SpecificHeat/model in input file.");
      }

    model = input( "Materials/"+material+"/SpecificHeat/model", "DIE!" );
  }

  void MaterialsParsing::turb_viscosity_model( const GetPot & input,
                                               const std::string & /*physics*/,
                                               const std::string & material,
                                               std::string & model )
  {
    if( !input.have_variable("Materials/"+material+"/Viscosity/turb_visc_model") )
      {
        libmesh_error_msg("Error: Could not find Materials/"+material+"/Viscosity/turb_visc_model in input file.");
      }

    model = input( "Materials/"+material+"/Viscosity/turb_visc_model", "DIE!" );
  }

  void MaterialsParsing::stress_strain_model( const GetPot & input,
                                              const std::string & /*physics*/,
                                              const std::string & material,
                                              std::string & model,
                                              std::string & strain_energy )
  {
    if( !input.have_variable("Materials/"+material+"/StressStrainLaw/model") )
      {
        libmesh_error_msg("Error: Could not find Materials/"+material+"/StressStrainLaw/model in input file!");
      }

    model = input("Materials/"+material+"/StressStrainLaw/model", "DIE!");
    strain_energy = input("Materials/"+material+"/StressStrainLaw/strain_energy", "none");
    // These options are deprecated
    if( model == std::string( "HookesLaw" ) )
      {
        std::string warning = "WARNING: Detected model HookesLaw. This is DEPRECATED!\n";
        warning += "         Please update to use hookes_law.\n";
        grins_warning(warning);

        model = std::string("hookes_law");
      }

    if( model == std::string( "MooneyRivlin" ) )
      {
        std::string warning = "WARNING: Detected model MooneyRivlin. This is DEPRECATED!\n";
        warning += "         Please update model to use incompressible_hyperelasticity and\n";
        warning += "         set strain_energy to mooney_rivlin.\n";
        grins_warning(warning);

        model = std::string("incompressible_hyperelasticity");
        strain_energy = std::string("mooney_rivlin");
      }
  }

  void MaterialsParsing::read_density( const std::string & core_physics_name,
                                       const GetPot & input,
                                       ParameterUser& params,
                                       libMesh::Real& rho )
  {
    bool have_material = MaterialsParsing::have_material(input,core_physics_name);

    std::string material = input("Physics/"+core_physics_name+"/material", "DIE!");

    // Error if both material/Density and rho are specified
    MaterialsParsing::duplicate_input_test(input,
                                           "Physics/"+core_physics_name+"/rho",
                                           "Materials/"+material+"/Density/value" );

    // It's deprecated to have nothing and default to 1.0
    if( !input.have_variable("Physics/"+core_physics_name+"/rho") &&
        !have_material )
      {
        // For some insane reason, we'd originally tied density specifically
        // to PhysicsNaming::incompressible_navier_stokes(), so we'll check for that first.
        if( input.have_variable("Physics/"+PhysicsNaming::incompressible_navier_stokes()+"/rho") )
          {
            std::string warning = "WARNING: neither Physics/"+core_physics_name+"/rho nor\n";
            warning += "         Physics/"+core_physics_name+"/material options were detected.\n";
            warning += "         But we found Physics/"+PhysicsNaming::incompressible_navier_stokes()+"/rho so we using that.\n";
            warning += "        This behavior is DEPRECATED.\n";
            warning += "         Please update and use Physics/"+core_physics_name+"/material.\n";
            grins_warning(warning);

            params.set_parameter
              (rho, input,
               "Physics/"+PhysicsNaming::incompressible_navier_stokes()+"/rho", 1.0 /*default*/);
          }
        // Otherwise, the insanity continued and we defaulted to 1.0
        else
          {
            std::string warning = "WARNING: neither Physics/"+core_physics_name+"/rho nor\n";
            warning += "         Physics/"+core_physics_name+"/material  nor\n";
            warning += "         Physics/"+PhysicsNaming::incompressible_navier_stokes()+"/rho options were detected.\n";
            warning += "         We are assuming a density value of 1.0. This is DEPRECATED.\n";
            warning += "         Please update and use Physics/"+core_physics_name+"/material.\n";
            grins_warning(warning);

            params.set_parameter
              (rho, input,
               "Physics/"+core_physics_name+"/rho", 1.0 /*default*/);
          }
      }
    // It's deprecated to use rho as the density input
    else if( input.have_variable("Physics/"+core_physics_name+"/rho") )
      {
        MaterialsParsing::dep_input_warning( "Physics/"+core_physics_name+"/rho",
                                             "Density/value" );

        params.set_parameter
          (rho, input,
           "Physics/"+core_physics_name+"/rho", 1.0 /*default*/);
      }

    // This is the preferred version
    else if( have_material )
      {
        if( !input.have_variable("Materials/"+material+"/Density/value") )
          libmesh_error_msg("ERROR: Could not find Materials/"+material+"/Density/value in input!");

        params.set_parameter
          (rho, input,
           "Materials/"+material+"/Density/value", 0.0 /*default*/);
      }
    // Wat?
    else
      libmesh_error();

    // Let's make sure we actually got a valid density value
    if( rho <= 0.0 )
      {
        libmesh_error_msg("ERROR: Detected non-positive value of density!");
      }
  }

  void MaterialsParsing::read_specific_heat( const std::string & core_physics_name,
                                             const GetPot & input,
                                             ParameterUser& params,
                                             libMesh::Real& cp )
  {
    std::string material = input("Physics/"+core_physics_name+"/material", "DIE!");

    // Error if both material/SpecificHeat and Cp are specified
    MaterialsParsing::duplicate_input_test(input,
                                           "Physics/"+core_physics_name+"/Cp",
                                           "Materials/"+material+"/SpecificHeat/value" );

    // It's deprecated to have nothing and default to 1.0
    if( !input.have_variable("Physics/"+core_physics_name+"/Cp") &&
        ( !input.have_variable("Physics/"+core_physics_name+"/material") ||
          !input.have_variable("Materials/"+material+"/SpecificHeat/value") ) )
      {
        // For some insane reason, we'd originally tied Cp specifically
        // to PhysicsNaming::heat_transfer(), so we'll check for that first.
        if( input.have_variable("Physics/"+PhysicsNaming::heat_transfer()+"/Cp") )
          {
            std::string warning = "WARNING: neither Physics/"+core_physics_name+"/Cp nor\n";
            warning += "         Physics/"+core_physics_name+"/material options were detected.\n";
            warning += "         But we found Physics/"+PhysicsNaming::heat_transfer()+"/Cp so we using that.\n";
            warning += "        This behavior is DEPRECATED.\n";
            warning += "         Please update and use Physics/"+core_physics_name+"/material.\n";
            grins_warning(warning);

            params.set_parameter
              (cp, input,
               "Physics/"+PhysicsNaming::heat_transfer()+"/Cp", 1.0 /*default*/);
          }
        // Otherwise, the insanity continued and we defaulted to 1.0
        else
          {
            std::string warning = "WARNING: neither Physics/"+core_physics_name+"/Cp nor\n";
            warning += "         Physics/"+core_physics_name+"/material  nor\n";
            warning += "         Physics/"+PhysicsNaming::heat_transfer()+"/Cp options were detected.\n";
            warning += "         We are assuming a specific heat value of 1.0. This is DEPRECATED.\n";
            warning += "         Please update and use Physics/"+core_physics_name+"/material.\n";
            grins_warning(warning);

            params.set_parameter
              (cp, input,
               "Physics/"+core_physics_name+"/Cp", 1.0 /*default*/);
          }
      }

    // It's deprecated to use Cp as the specific heat input
    if( input.have_variable("Physics/"+core_physics_name+"/Cp") )
      {
        MaterialsParsing::dep_input_warning( "Physics/"+core_physics_name+"/Cp",
                                             "SpecificHeat/value" );

        params.set_parameter
          (cp, input,
           "Physics/"+core_physics_name+"/Cp", 1.0 /*default*/);
      }

    // This is the preferred version
    if( input.have_variable("Physics/"+core_physics_name+"/material") &&
        input.have_variable("Materials/"+material+"/SpecificHeat/value") )
      {
        // This function only supports reading a constant
        if( input("Materials/"+material+"/SpecificHeat/model", "DIE!") !=
            std::string("constant") )
          libmesh_error_msg("ERROR: Only constant SpecificHeat model supported!");

        params.set_parameter
          (cp, input,
           "Materials/"+material+"/SpecificHeat/value", 0.0 /*default*/);
      }

    // Let's make sure we actually got a valid density value
    if( cp <= 0.0 )
      {
        libmesh_error_msg("ERROR: Detected non-positive value of cp!");
      }
  }

  void MaterialsParsing::read_property( const GetPot & input,
                                        const std::string & old_option,
                                        const std::string & property,
                                        const std::string & core_physics,
                                        ParameterUser& param_user,
                                        libMesh::Real& value )
  {
    std::string material = MaterialsParsing::material_name(input,core_physics);

    // Can't specify both old_option and property
    MaterialsParsing::duplicate_input_test(input,
                                           old_option,
                                           "Materials/"+material+"/"+property+"/value" );

    // Deprecated
    if( input.have_variable(old_option) )
      {
        MaterialsParsing::dep_input_warning( old_option,property+"/value" );

        param_user.set_parameter(value, input, old_option, 0.0 /*default*/);
      }
    // Preferred
    else if( input.have_variable("Materials/"+material+"/"+property+"/value" ) )
      {
        param_user.set_parameter
          (value, input, "Materials/"+material+"/"+property+"/value", 0.0 /*default*/);
      }
    // If nothing was set, that's an error
    else
      {
        libmesh_error_msg("ERROR: No valid input found for "+property+"!");
      }

    // Make sure value is positive
    if( value <= 0.0 )
      {
        libmesh_error_msg("ERROR: Detected non-positive "+property+"!");
      }
  }

  void MaterialsParsing::parse_chemical_species( const GetPot & input,
                                                 const std::string & material,
                                                 std::vector<std::string>& species_names )
  {
    // Clear out anything the user might've put in there.
    species_names.clear();

    MaterialsParsing::duplicate_input_test( input,
                                            "Physics/Chemistry/species",
                                            "Materials/"+material+"/GasMixture/species");

    std::string species_input;
    if( input.have_variable("Physics/Chemistry/species") )
      {
        MaterialsParsing::dep_input_warning("Physics/Chemistry/species","GasMixture/species");
        species_input = "Physics/Chemistry/species";
      }
    else if( input.have_variable("Materials/"+material+"/GasMixture/species") )
      {
        species_input = "Materials/"+material+"/GasMixture/species";
      }
    else
      {
        libmesh_error_msg("ERROR: Valid input for species not found!");
      }

    // Read variable naming info
    unsigned int n_species = input.vector_variable_size(species_input);

    species_names.reserve(n_species);
    for( unsigned int i = 0; i < n_species; i++ )
      {
        /*! \todo Make this prefix string an input option */

        species_names.push_back( input( species_input, "DIE!", i ) );
      }
  }

  void MaterialsParsing::parse_species_varnames( const GetPot & input,
                                                 const std::string & material,
                                                 const std::string & prefix,
                                                 std::vector<std::string>& species_varnames )
  {
    std::vector<std::string> species_names;
    MaterialsParsing::parse_chemical_species(input,material,species_names);
    unsigned int n_species = species_names.size();
    species_varnames.reserve(n_species);

    for( unsigned int i = 0; i < n_species; i++ )
      {
        std::string var_name = prefix+species_names[i];
        species_varnames.push_back(var_name);
      }
  }

  libMesh::Real MaterialsParsing::parse_lewis_number( const GetPot & input,
                                                      const std::string & material )
  {

    // Can't specify both old_option and property
    MaterialsParsing::duplicate_input_test(input,
                                           "Physics/Antioch/Le",
                                           "Materials/"+material+"/LewisNumber/value");

    libMesh::Real Le = 0.0;
    if( input.have_variable("Physics/Antioch/Le") )
      {
        MaterialsParsing::dep_input_warning("Physics/Antioch/Le", "LewisNumber/value" );
        Le = input("Physics/Antioch/Le", 0.0);
      }
    else if( input.have_variable("Materials/"+material+"/LewisNumber/value") )
      {
        Le = input("Materials/"+material+"/LewisNumber/value", 0.0);
      }
    else
      {
        libmesh_error_msg("ERROR: Could not find value input for LewisNumber!");
      }

    return Le;
  }

  std::string MaterialsParsing::parse_chemical_kinetics_datafile_name( const GetPot & input,
                                                                       const std::string & material )
  {
    // Can't specify both old_option and property
    MaterialsParsing::duplicate_input_test(input,
                                           "Physics/Chemistry/chem_file",
                                           "Materials/"+material+"/GasMixture/kinetics_data");

    std::string filename;
    if( input.have_variable("Physics/Chemistry/chem_file") )
      {
        MaterialsParsing::dep_input_warning("Physics/Chemistry/chem_file",
                                            "GasMixture/kinetics_data" );

        filename = input("Physics/Chemistry/chem_file","DIE!");
      }
    else if( input.have_variable("Materials/"+material+"/GasMixture/kinetics_data") )
      {
        filename = input("Materials/"+material+"/GasMixture/kinetics_data", "DIE!");
      }
    else
      {
        libmesh_error_msg("ERROR: Could not find valid input for kinetics_data!");
      }

    return filename;
  }

  void MaterialsParsing::dep_input_warning( const std::string & old_option,
                                            const std::string & property )
  {
    std::string warning = "WARNING: Input option "+old_option+" is DEPRECATED!\n";
    warning += "         Please update to use Materials/MATERIAL_NAME/"+property+"\n";
    grins_warning(warning);
  }

  void MaterialsParsing::duplicate_input_test( const GetPot & input,
                                               const std::string & option1,
                                               const std::string & option2 )
  {
    if( input.have_variable(option1) &&
        input.have_variable(option2) )
      {
        libmesh_error_msg("ERROR: Can't specify both "+option1+" and "+option2+"!");
      }
  }

} // end namespace GRINS
