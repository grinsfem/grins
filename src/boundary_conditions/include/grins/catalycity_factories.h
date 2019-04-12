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

#ifndef GRINS_CATALYCITY_FACTORIES_H
#define GRINS_CATALYCITY_FACTORIES_H

// C+
#include <limits>

// GRINS
#include "grins/catalycity_factory_abstract.h"
#include "grins/constant_catalycity.h"
#include "grins/arrhenius_catalycity.h"
#include "grins/power_law_catalycity.h"

namespace GRINS
{
  class ConstantCatalycityFactory : public CatalycityFactoryAbstract
  {
  public:

    ConstantCatalycityFactory( const std::string& physics_name )
      : CatalycityFactoryAbstract(physics_name)
    {}

    ~ConstantCatalycityFactory(){};

  protected:

    virtual std::unique_ptr<CatalycityBase>
    build_catalycity( const GetPot& input, const std::string& section )
    {
      std::string param_base = section+"/ConstantCatalycity/";

      std::string gamma_str = param_base+"gamma";
      if( !input.have_variable(gamma_str) )
        libmesh_error_msg("ERROR: Could not find input "+gamma_str+" for ConstantCatalycity!\n");

      libMesh::Real gamma = input(gamma_str, std::numeric_limits<libMesh::Real>::max());
      std::unique_ptr<CatalycityBase> catalycity( new ConstantCatalycity( gamma ) );
      catalycity->set_parameters(input,param_base);
      return catalycity;
    }

  };


  class ArrheniusCatalycityFactory : public CatalycityFactoryAbstract
  {
  public:

    ArrheniusCatalycityFactory( const std::string& physics_name )
      : CatalycityFactoryAbstract(physics_name)
    {}

    ~ArrheniusCatalycityFactory(){};

  protected:

    virtual std::unique_ptr<CatalycityBase>
    build_catalycity( const GetPot& input, const std::string& section )
    {
      std::string param_base = section+"/ArrheniusCatalycity/";

      std::string gamma_str = param_base+"gamma0";
      if( !input.have_variable(gamma_str) )
        libmesh_error_msg("ERROR: Could not find input "+gamma_str+" for ArrheniusCatalycity!\n");

      std::string Ta_str = param_base+"Ta";
      if( !input.have_variable(Ta_str) )
        libmesh_error_msg("ERROR: Could not find input "+Ta_str+" for ArrheniusCatalycity!\n");

      libMesh::Real gamma = input(gamma_str, std::numeric_limits<libMesh::Real>::max());
      libMesh::Real Ta = input(Ta_str, std::numeric_limits<libMesh::Real>::max());

      std::unique_ptr<CatalycityBase> catalycity( new ArrheniusCatalycity( gamma, Ta ) );
      catalycity->set_parameters(input,param_base);
      return catalycity;
    }

  };

  class PowerLawCatalycityFactory : public CatalycityFactoryAbstract
  {
  public:

    PowerLawCatalycityFactory( const std::string& physics_name )
      : CatalycityFactoryAbstract(physics_name)
    {}

    ~PowerLawCatalycityFactory(){};

  protected:

    virtual std::unique_ptr<CatalycityBase>
    build_catalycity( const GetPot& input, const std::string& section )
    {
      std::string param_base = section+"/PowerLawCatalycity/";

      std::string gamma_str = param_base+"gamma0";
      if( !input.have_variable(gamma_str) )
        libmesh_error_msg("ERROR: Could not find input "+gamma_str+" for ArrheniusCatalycity!\n");

      std::string Tref_str = param_base+"Tref";
      if( !input.have_variable(Tref_str) )
        libmesh_error_msg("ERROR: Could not find input "+Tref_str+" for PowerLawCatalycity!\n");

      std::string alpha_str = param_base+"alpha";
      if( !input.have_variable(alpha_str) )
        libmesh_error_msg("ERROR: Could not find input "+alpha_str+" for PowerLawCatalycity!\n");

      libMesh::Real gamma = input(gamma_str, std::numeric_limits<libMesh::Real>::max());
      libMesh::Real Tref = input(Tref_str, std::numeric_limits<libMesh::Real>::max());
      libMesh::Real alpha = input(alpha_str, std::numeric_limits<libMesh::Real>::max());

      std::unique_ptr<CatalycityBase> catalycity( new PowerLawCatalycity( gamma, Tref, alpha ) );
      catalycity->set_parameters(input,param_base);
      return catalycity;
    }

  };

} // end namespace GRINS

#endif // GRINS_CATALYCITY_FACTORIES_H
