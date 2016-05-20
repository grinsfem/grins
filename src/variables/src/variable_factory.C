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

// This class
#include "grins/variable_factory.h"

namespace GRINS
{
  libMesh::UniquePtr<FEVariablesBase> VariableFactoryAbstract::create()
  {
    // Make sure all necessary state has been setup
    this->check_create_state();

    libMesh::UniquePtr<FEVariablesBase> func =
      this->build_fe_var( *_var_names, *_var_indices );

    // Reset state for error checking during next construction
    this->reset_create_state();

    return func;
  }

  std::vector<std::string> VariableFactoryAbstract::build_var_names( const std::string& name )
  {
    if( !_input )
      libmesh_error_msg("ERROR: Must call set_getpot() before calling VariableFactoryAbstract::build_var_names!");

    if( _var_section == std::string("DIE!") )
      libmesh_error_msg("ERROR: Must call set_var_section() before calling VariableFactoryAbstract::build_var_names!");

    VariableFactoryAbstract& factory = get_factory_subclass<VariableFactoryAbstract>(name);

    std::vector<std::string> var_names;
    var_names = factory.parse_var_names( *_input, _var_section );

    // Reset _input to NULL for error checking
    _input = NULL;
    _var_section = std::string("DIE!");

    return var_names;
  }

  void VariableFactoryAbstract::check_create_state() const
  {
    if( !this->_var_names )
      libmesh_error_msg("ERROR: must call set_var_names() before building FEVariablesBase!");

    if( !this->_var_indices )
      libmesh_error_msg("ERROR: must call set_var_indices() before building FEVariablesBase!");
  }

  void VariableFactoryAbstract::reset_create_state()
  {
    _var_names = NULL;
    _var_indices = NULL;
  }


  // Full specialization for the Factory<FEVariablesBase>
  template<>
  std::map<std::string, FactoryAbstract<FEVariablesBase>*>&
  FactoryAbstract<FEVariablesBase>::factory_map()
  {
    static std::map<std::string, FactoryAbstract<FEVariablesBase>*> _map;
    return _map;
  }

  // Definition of static members
  template<>
  const GetPot* FactoryWithGetPot<FEVariablesBase>::_input = NULL;
  const std::vector<std::string>* VariableFactoryAbstract::_var_names = NULL;
  const std::vector<VariableIndex>* VariableFactoryAbstract::_var_indices = NULL;
  std::string VariableFactoryAbstract::_var_section = std::string("DIE!");

} // end namespace GRINS
