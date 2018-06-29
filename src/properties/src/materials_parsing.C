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
                                              const std::string & physics,
                                              std::string & model )
  {
    std::string material = MaterialsParsing::material_name( input, physics );
    std::string spec_heat_option("Materials/"+material+"/SpecificHeat/model");

    MaterialsParsing::check_for_input_option(input,spec_heat_option);

    model = input( spec_heat_option, "DIE!" );
  }

  void MaterialsParsing::turb_viscosity_model( const GetPot & input,
                                               const std::string & physics,
                                               std::string & model )
  {
    std::string material = MaterialsParsing::material_name( input, physics );
    std::string visc_option("Materials/"+material+"/Viscosity/turb_visc_model");

    MaterialsParsing::check_for_input_option(input,visc_option);

    model = input( visc_option, "DIE!" );
  }

  void MaterialsParsing::stress_strain_model( const GetPot & input,
                                              const std::string & physics,
                                              std::string & model,
                                              std::string & strain_energy )
  {
    std::string material = MaterialsParsing::material_name( input, physics );
    std::string ss_option("Materials/"+material+"/StressStrainLaw/model");

    MaterialsParsing::check_for_input_option(input,ss_option);

    model = input(ss_option, "DIE!");

    strain_energy = input("Materials/"+material+"/StressStrainLaw/strain_energy", "none");
  }

  void MaterialsParsing::thermochemistry_lib( const GetPot & input,
                                                const std::string & physics,
                                                std::string & thermochem_lib )
  {
    std::string material = MaterialsParsing::material_name( input, physics );

    std::string thermochem_lib_option =
      MaterialsParsing::thermochem_lib_input_string(material);

    MaterialsParsing::check_for_input_option(input,thermochem_lib_option);

    thermochem_lib = input(thermochem_lib_option, "DIE!");

    // Make sure we have a valid model
    if( thermochem_lib != std::string("antioch") &&
        thermochem_lib != std::string("cantera") )
      {
        std::string error = "ERROR! Invalid thermochemistry library value "+thermochem_lib+"!\n";
        error += "       Valid selections are: antioch\n";
        error += "                             cantera\n";

        libmesh_error_msg(error);
      }
  }

  void MaterialsParsing::antioch_models( const GetPot& input,
                                         const std::string& physics,
                                         std::string& transport_model,
                                         std::string& thermo_model,
                                         std::string& viscosity_model,
                                         std::string& conductivity_model,
                                         std::string& diffusivity_model )
  {
    // Newer, preferred version
    std::string material = MaterialsParsing::material_name( input, physics );

    MaterialsParsing::check_for_input_option(input,"Materials/"+material+"/GasMixture/Antioch/transport_model");
    transport_model = input("Materials/"+material+"/GasMixture/Antioch/transport_model", "DIE!");

    // Now parse the remaining models
    MaterialsParsing::check_for_input_option(input,"Materials/"+material+"/GasMixture/Antioch/thermo_model");
    MaterialsParsing::check_for_input_option(input,"Materials/"+material+"/GasMixture/Antioch/viscosity_model");
    MaterialsParsing::check_for_input_option(input,"Materials/"+material+"/GasMixture/Antioch/thermal_conductivity_model");
    MaterialsParsing::check_for_input_option(input,"Materials/"+material+"/GasMixture/Antioch/mass_diffusivity_model");

    thermo_model = input( "Materials/"+material+"/GasMixture/Antioch/thermo_model", "DIE!");
    viscosity_model = input( "Materials/"+material+"/GasMixture/Antioch/viscosity_model", "DIE!");
    conductivity_model = input( "Materials/"+material+"/GasMixture/Antioch/thermal_conductivity_model", "DIE!");
    diffusivity_model = input( "Materials/"+material+"/GasMixture/Antioch/mass_diffusivity_model", "DIE!");
  }

  void MaterialsParsing::read_property( const GetPot & input,
                                        const std::string & property,
                                        const std::string & core_physics,
                                        ParameterUser& param_user,
                                        libMesh::Real& value )
  {
    std::string material = MaterialsParsing::material_name(input,core_physics);

    std::string option("Materials/"+material+"/"+property+"/value");
    MaterialsParsing::check_for_input_option(input,option);

    param_user.set_parameter(value, input, option, 0.0 /*default*/);

    if( value <= 0.0 )
      libmesh_error_msg("ERROR: Detected non-positive "+property+"!");
  }

  void MaterialsParsing::parse_chemical_species( const GetPot & input,
                                                 const std::string & material,
                                                 std::vector<std::string>& species_names )
  {
    // Clear out anything the user might've put in there.
    species_names.clear();

    std::string option("Materials/"+material+"/GasMixture/species");
    MaterialsParsing::check_for_input_option(input,option);

    // Read variable naming info
    unsigned int n_species = input.vector_variable_size(option);

    species_names.reserve(n_species);
    for( unsigned int i = 0; i < n_species; i++ )
      species_names.push_back( input( option, "DIE!", i ) );
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

    std::string option("Materials/"+material+"/LewisNumber/value");
    MaterialsParsing::check_for_input_option(input,option);

    libMesh::Real Le = input(option, 0.0);

    return Le;
  }

  std::string MaterialsParsing::parse_chemical_kinetics_datafile_name( const GetPot & input,
                                                                       const std::string & material )
  {
    std::string option("Materials/"+material+"/GasMixture/kinetics_data");
    MaterialsParsing::check_for_input_option(input,option);

    std::string warning = "WARNING: option "+option+"is DEPRECATED!\n";
    warning += "         kinetics_data moved to thermochemistry section\n";
    warning += "         and renamed chemical_data!";
    grins_warning_once(warning);

    std::string filename = input(option, "DIE!");

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
