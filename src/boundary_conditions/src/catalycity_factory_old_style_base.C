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

// This class
#include "grins/catalycity_factory_old_style_base.h"

namespace GRINS
{
  std::unique_ptr<CatalycityBase> CatalycityFactoryOldStyleBase::build_catalycity( const GetPot& input,
                                                                                   const std::string& section )
  {
    // State of _reactant_str, _bc_id_str verified in check_state() call in create()
    return this->build_catalycity_old_style( input, section, _reactant_str, _bc_id_str );
  }

  void CatalycityFactoryOldStyleBase::check_state() const
  {
    CatalycityFactoryAbstract::check_state();

    if( _reactant_str == std::string("DIE!") )
      libmesh_error_msg("ERROR: must call set_reactant() before building Catalycity!");

    if( _bc_id_str == std::string("DIE!") )
      libmesh_error_msg("ERROR: must call set_bc_id() before building Catalycity!");
  }

  void CatalycityFactoryOldStyleBase::reset_state()
  {
    CatalycityFactoryAbstract::reset_state();

    _reactant_str = std::string("DIE!");
    _bc_id_str = std::string("DIE!");
  }

  // Definition of static members
  std::string CatalycityFactoryOldStyleBase::_reactant_str = std::string("DIE!");
  std::string CatalycityFactoryOldStyleBase::_bc_id_str = std::string("DIE!");

} // end namespace GRINS
