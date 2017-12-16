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

// This class
#include "grins/source_term_base.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  SourceTermBase::SourceTermBase( const std::string& physics_name, const GetPot& input )
    : Physics(physics_name,input)
  {
    this->parse_var_info(input);

    return;
  }

  SourceTermBase::~SourceTermBase()
  {
    return;
  }

  void SourceTermBase::parse_var_info( const GetPot& input )
  {
    if( !input.have_variable("Physics/"+this->_physics_name+"/Variables/names") )
      {
        libMesh::err << "Error: Must have at least one variable for source function." << std::endl
                     << "       Ensure that Physics/"+this->_physics_name+"/Variables/names is set." << std::endl;
        libmesh_error();
      }

    unsigned int n_vars = input.vector_variable_size("Physics/"+this->_physics_name+"/Variables/names");

    // Make sure we have consisent number of FE types and FE orders
    /*! \todo In the future, after refactoring Variable parsing, we should be
      able get the FE type and order information from there. */
    unsigned int n_fe_types = input.vector_variable_size("Physics/"+this->_physics_name+"/Variables/FE_types");

    unsigned int n_fe_orders = input.vector_variable_size("Physics/"+this->_physics_name+"/Variables/FE_orders");

    if( n_fe_types != n_vars )
      {
        libMesh::err << "Error: Must have matching number of variable names and FE types." << std::endl
                     << "       Found " << n_fe_types << " FE types and " << n_vars << " variables." << std::endl
                     << "       Ensure Physics/"+this->_physics_name+"/Variables/FE_types is consistent." << std::endl;
        libmesh_error();
      }

    if( n_fe_orders != n_vars )
      {
        libMesh::err << "Error: Must have matching number of variable names and FE orders." << std::endl
                     << "       Found " << n_fe_orders << " FE orders and " << n_vars << " variables." << std::endl
                     << "       Ensure Physics/"+this->_physics_name+"/Variables/FE_orders is consistent." << std::endl;
        libmesh_error();
      }

    _var_names.reserve(n_vars);
    _var_FE.reserve(n_vars);
    _var_order.reserve(n_vars);
    for( unsigned int v = 0; v < n_vars; v++ )
      {
        _var_names.push_back( input("Physics/"+this->_physics_name+"/Variables/names", "DIE!", v) );
        _var_FE.push_back( libMesh::Utility::string_to_enum<GRINSEnums::FEFamily>(input("Physics/"+this->_physics_name+"/Variables/FE_types", "DIE!", v)) );
        _var_order.push_back( libMesh::Utility::string_to_enum<GRINSEnums::Order>(input("Physics/"+this->_physics_name+"/Variables/FE_orders", "DIE!", v)) );
      }

    return;
  }

  void SourceTermBase::init_variables( libMesh::FEMSystem* system )
  {
    // We'd better have at least 1 variable read from input
    libmesh_assert( !_var_names.empty() );

    _vars.resize( _var_names.size() );

    for( unsigned int var = 0; var < _vars.size(); var++ )
      {
        _vars[var] = system->add_variable( _var_names[var], _var_order[var], _var_FE[var] );
      }

    return;
  }
} // end namespace GRINS
