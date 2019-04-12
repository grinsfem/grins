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

#ifndef GRINS_CATALYCITY_FACTORIES_OLD_STYLE_H
#define GRINS_CATALYCITY_FACTORIES_OLD_STYLE_H

// C+
#include <limits>

// GRINS
#include "grins/catalycity_factory_old_style_base.h"
#include "grins/constant_catalycity.h"
#include "grins/arrhenius_catalycity.h"
#include "grins/power_law_catalycity.h"

namespace GRINS
{
  class ConstantCatalycityFactoryOldStyle : public CatalycityFactoryOldStyleBase
  {
  public:

    ConstantCatalycityFactoryOldStyle( const std::string& physics_name )
      : CatalycityFactoryOldStyleBase(physics_name)
    {}

    ~ConstantCatalycityFactoryOldStyle(){};

  protected:

    virtual std::unique_ptr<CatalycityBase> build_catalycity_old_style( const GetPot& input,
                                                                        const std::string& section,
                                                                        const std::string& reactant_str,
                                                                        const std::string& bc_id_string )
    {
      std::string gamma_str = section+"/gamma_"+reactant_str+"_"+bc_id_string;
      if( !input.have_variable(gamma_str) )
        libmesh_error_msg("ERROR: Could not find input "+gamma_str+" for ConstantCatalycity!\n");

      libMesh::Real gamma = input(gamma_str, std::numeric_limits<libMesh::Real>::max());
      return std::unique_ptr<CatalycityBase>( new ConstantCatalycity( gamma ) );
    }
  };

  class ArrheniusCatalycityFactoryOldStyle : public CatalycityFactoryOldStyleBase
  {
  public:

    ArrheniusCatalycityFactoryOldStyle( const std::string& physics_name )
      : CatalycityFactoryOldStyleBase(physics_name)
    {}

    ~ArrheniusCatalycityFactoryOldStyle(){};

  protected:

    virtual std::unique_ptr<CatalycityBase> build_catalycity_old_style( const GetPot& input,
                                                                        const std::string& section,
                                                                        const std::string& reactant_str,
                                                                        const std::string& bc_id_string )
    {
      std::string gamma_str = section+"/gamma0_"+reactant_str+"_"+bc_id_string;
      if( !input.have_variable(gamma_str) )
        libmesh_error_msg("ERROR: Could not find input "+gamma_str+" for ArrheniusCatalycity!\n");

      std::string Ta_str = section+"/Ta_"+reactant_str+"_"+bc_id_string;
      if( !input.have_variable(Ta_str) )
        libmesh_error_msg("ERROR: Could not find input "+Ta_str+" for ArrheniusCatalycity!\n");

      libMesh::Real gamma = input(gamma_str, std::numeric_limits<libMesh::Real>::max());
      libMesh::Real Ta = input(Ta_str, std::numeric_limits<libMesh::Real>::max());

      return std::unique_ptr<CatalycityBase>( new ArrheniusCatalycity( gamma, Ta ) );
    }
  };

  class PowerLawCatalycityFactoryOldStyle : public CatalycityFactoryOldStyleBase
  {
  public:

    PowerLawCatalycityFactoryOldStyle( const std::string& physics_name )
      : CatalycityFactoryOldStyleBase(physics_name)
    {}

    ~PowerLawCatalycityFactoryOldStyle(){};

  protected:

    virtual std::unique_ptr<CatalycityBase> build_catalycity_old_style( const GetPot& input,
                                                                        const std::string& section,
                                                                        const std::string& reactant_str,
                                                                        const std::string& bc_id_string )
    {
      std::string gamma_str = section+"/gamma0_"+reactant_str+"_"+bc_id_string;
      if( !input.have_variable(gamma_str) )
        libmesh_error_msg("ERROR: Could not find input "+gamma_str+" for PowerLawCatalycity!\n");

      std::string Tref_str = section+"/Tref_"+reactant_str+"_"+bc_id_string;
      if( !input.have_variable(Tref_str) )
        libmesh_error_msg("ERROR: Could not find input "+Tref_str+" for PowerLawCatalycity!\n");

      std::string alpha_str = section+"/alpha_"+reactant_str+"_"+bc_id_string;
      if( !input.have_variable(alpha_str) )
        libmesh_error_msg("ERROR: Could not find input "+alpha_str+" for PowerLawCatalycity!\n");

      libMesh::Real gamma = input(gamma_str, std::numeric_limits<libMesh::Real>::max());
      libMesh::Real Tref = input(Tref_str, std::numeric_limits<libMesh::Real>::max());
      libMesh::Real alpha = input(alpha_str, std::numeric_limits<libMesh::Real>::max());

      return std::unique_ptr<CatalycityBase>( new PowerLawCatalycity( gamma, Tref, alpha ) );
    }
  };

} // end namespace GRINS

#endif // GRINS_CATALYCITY_FACTORIES_OLD_STYLE_H
