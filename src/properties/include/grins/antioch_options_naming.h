//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_ANTIOCH_OPTIONS_NAMING_H
#define GRINS_ANTIOCH_OPTIONS_NAMING_H

namespace GRINS
{
  class AntiochOptions
  {
  public:
    AntiochOptions(){}
    ~AntiochOptions(){}

    static std::string stat_mech_thermo_model()
    { return "stat_mech";}

    static std::string ideal_gas_thermo_model()
    { return "ideal_gas";}

    static std::string cea_nasa_model()
    { return "cea";}

    static std::string nasa7_nasa_model()
    { return "nasa7";}

    static std::string nasa9_nasa_model()
    { return "nasa9";}

    static std::string constant_transport_model()
    { return "constant";}

    static std::string mix_avged_transport_model()
    { return "mixture_averaged";}

    static std::string constant_conductivity_model()
    { return "constant";}

    static std::string constant_prandtl_conductivity_model()
    { return "constant_prandtl";}

    static std::string eucken_conductivity_model()
    { return "eucken";}

    static std::string kinetic_theory_conductivity_model()
    { return "kinetics_theory";}

    static std::string constant_viscosity_model()
    { return "constant";}

    static std::string sutherland_viscosity_model()
    { return "sutherland";}

    static std::string blottner_viscosity_model()
    { return "blottner";}

    static std::string kinetic_theory_viscosity_model()
    { return "kinetics_theory";}

    static std::string constant_lewis_diffusivity_model()
    { return "constant_lewis";}

    static std::string kinetic_theory_diffusivity_model()
    { return "kinetics_theory";}

  };
} // end namespace GRINS
#endif // GRINS_ANTIOCH_OPTIONS_NAMING_H
