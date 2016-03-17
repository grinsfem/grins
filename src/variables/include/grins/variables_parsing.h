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

#ifndef GRINS_VARIABLES_PARSING_H
#define GRINS_VARIABLES_PARSING_H

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

  };

} // end namespace GRINS

#endif // GRINS_VARIABLES_PARSING_H
