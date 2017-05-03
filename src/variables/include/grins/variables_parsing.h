//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_VARIABLES_PARSING_H
#define GRINS_VARIABLES_PARSING_H

// C++
#include <string>

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  class VariablesParsing
  {
  public:

    static std::string displacement_section()
    { return "Displacement"; }

    static std::string pressure_section()
    { return "Pressure"; }

    static std::string temperature_section()
    { return "Temperature"; }

    static std::string thermo_pressure_section()
    { return "ThermoPressure"; }

    static std::string turbulence_section()
    { return "TurbulentViscosity"; }

    static std::string species_mass_fractions_section()
    { return "SpeciesMassFractions"; }

    static std::string velocity_section()
    { return "Velocity"; }

    static std::string single_var_section()
    { return "SingleVariable"; }

    static std::string scalar_var_section()
    { return "ScalarVariable"; }

    //! Helper function to encapsualte the overall [Variables] section name.
    static std::string variables_section()
    { return "Variables"; }

    //! Helper function to encapsulate the names variable in the input file
    /*! This is the full variable name to be passed to GetPot to read in
      user-supplied names for variables. */
    static std::string varnames_input_name( const std::string& subsection )
    { return VariablesParsing::variables_section()+"/"+subsection+"/names"; }

    //! Helper function to encaplusate fe_family input variable
    static std::string fe_family_input_name( const std::string& subsection )
    { return VariablesParsing::variables_section()+"/"+subsection+"/fe_family"; }

    //! Helper function to encaplusate order input variable
    static std::string order_input_name( const std::string& subsection )
    { return VariablesParsing::variables_section()+"/"+subsection+"/order"; }

    enum SECTION_TYPE { PHYSICS = 0 };

    static std::string single_variable_name( const GetPot& input, const std::string& subsection_name,
                                             const SECTION_TYPE section_type )
    { return VariablesParsing::section_parse_var_name(input,
                                                      subsection_name,
                                                      "var_name",
                                                      VariablesParsing::single_var_section(),
                                                      section_type); }

    static std::string physics_scalar_variable_name( const GetPot& input, const std::string& physics_name )
    { return VariablesParsing::section_parse_var_name(input,
                                                      physics_name,
                                                      "var_name",
                                                      VariablesParsing::scalar_var_section(),
                                                      PHYSICS); }

    static std::string physics_velocity_variable_name( const GetPot& input, const std::string& physics_name )
    { return VariablesParsing::section_parse_var_name(input,
                                                      physics_name,
                                                      "velocity_var_name",
                                                      VariablesParsing::velocity_section(),
                                                      PHYSICS); }

    static std::string physics_temp_variable_name( const GetPot& input, const std::string& physics_name )
    { return VariablesParsing::section_parse_var_name(input,
                                                      physics_name,
                                                      "temperature_var_name",
                                                      VariablesParsing::temperature_section(),
                                                      PHYSICS); }

    static std::string physics_press_variable_name( const GetPot& input, const std::string& physics_name )
    { return VariablesParsing::section_parse_var_name(input,
                                                      physics_name,
                                                      "pressure_var_name",
                                                      VariablesParsing::pressure_section(),
                                                      PHYSICS); }

    static std::string physics_thermo_press_variable_name( const GetPot& input, const std::string& physics_name )
    { return VariablesParsing::section_parse_var_name(input,
                                                      physics_name,
                                                      "thermo_pressure_var_name",
                                                      VariablesParsing::thermo_pressure_section(),
                                                      PHYSICS); }

    static std::string physics_turb_variable_name( const GetPot& input, const std::string& physics_name )
    { return VariablesParsing::section_parse_var_name(input,
                                                      physics_name,
                                                      "turbulence_var_name",
                                                      VariablesParsing::turbulence_section(),
                                                      PHYSICS); }

    static std::string physics_disp_variable_name( const GetPot& input, const std::string& physics_name )
    { return VariablesParsing::section_parse_var_name(input,
                                                      physics_name,
                                                      "displacement_var_name",
                                                      VariablesParsing::displacement_section(),
                                                      PHYSICS); }

    static std::string physics_species_mass_frac_variable_name( const GetPot& input,
                                                                const std::string& physics_name )
    { return VariablesParsing::section_parse_var_name(input,
                                                      physics_name,
                                                      "species_mass_fracs_var_name",
                                                      VariablesParsing::species_mass_fractions_section(),
                                                      PHYSICS); }



  private:

    static std::string section_parse_var_name( const GetPot & input,
                                               const std::string & input_subsection_name,
                                               const std::string & input_var_name,
                                               const std::string & default_var_name,
                                               const SECTION_TYPE section_type );
  };

} // end namespace GRINS

#endif // GRINS_VARIABLES_PARSING_H
