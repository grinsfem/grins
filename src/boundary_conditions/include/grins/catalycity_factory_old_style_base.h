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

#ifndef GRINS_CATALYCITY_FACTORY_OLD_STYLE_BASE_H
#define GRINS_CATALYCITY_FACTORY_OLD_STYLE_BASE_H

// GRINS
#include "grins/catalycity_factory_abstract.h"

namespace GRINS
{
  class CatalycityFactoryOldStyleBase : public CatalycityFactoryAbstract
  {
  public:

    CatalycityFactoryOldStyleBase( const std::string& physics_name )
      : CatalycityFactoryAbstract(physics_name)
    {}

    ~CatalycityFactoryOldStyleBase(){};

    static void set_reactant( const std::string& reactant )
    { _reactant_str = reactant; }

    static void set_bc_id( const std::string& bc_id )
    { _bc_id_str = bc_id; }

  protected:

    virtual std::unique_ptr<CatalycityBase> build_catalycity( const GetPot& input,
                                                              const std::string& section );

    virtual std::unique_ptr<CatalycityBase> build_catalycity_old_style( const GetPot& input,
                                                                        const std::string& section,
                                                                        const std::string& reactant_str,
                                                                        const std::string& bc_id_string ) =0;
    virtual void check_state() const;

    virtual void reset_state();

    static std::string _reactant_str;

    static std::string _bc_id_str;

  };

} // end namespace GRINS

#endif // GRINS_CATALYCITY_FACTORY_OLD_STYLE_BASE_H
