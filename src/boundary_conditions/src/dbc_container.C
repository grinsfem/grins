//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "grins/dbc_container.h"

namespace GRINS
{

  DBCContainer::DBCContainer()
    : _var_names( std::vector<VariableName>() ),
      _bc_ids( std::set<BoundaryID>() ),
      _func( std::tr1::shared_ptr<libMesh::FunctionBase<libMesh::Number> >() )
  {
    return;
  }

  DBCContainer::~DBCContainer()
  {
    return;
  }

  void DBCContainer::add_var_name( const VariableName& var )
  {
    _var_names.push_back( var );
    return;
  }

  void DBCContainer::add_bc_id( const BoundaryID bc_id )
  {
    _bc_ids.insert( bc_id );
    return;
  }

  void DBCContainer::set_func( std::tr1::shared_ptr<libMesh::FunctionBase<libMesh::Number> > func )
  {
    _func = func;
    return;
  }

  std::vector<VariableName> DBCContainer::get_var_names() const
  {
    return _var_names;
  }

  std::set<BoundaryID> DBCContainer::get_bc_ids() const
  {
    return _bc_ids;
  }

  std::tr1::shared_ptr<libMesh::FunctionBase<libMesh::Number> > DBCContainer::get_func() const
  {
    return _func;
  }

} // namespace GRINS
