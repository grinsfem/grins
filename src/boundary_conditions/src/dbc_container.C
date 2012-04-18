//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "dbc_container.h"

GRINS::DBCContainer::DBCContainer()
  : _var_names( std::vector<GRINS::VariableName>() ),
    _bc_ids( std::set<GRINS::BoundaryID>() ),
    _func( std::tr1::shared_ptr<libMesh::FunctionBase<Number> >() )
{
  return;
}

GRINS::DBCContainer::~DBCContainer()
{
  return;
}

void GRINS::DBCContainer::set_var_names( const std::vector<GRINS::VariableName>& var_names )
{
  _var_names = var_names;
  return;
}

void GRINS::DBCContainer::set_bc_ids( const std::set<GRINS::BoundaryID>& bc_ids )
{
  _bc_ids = bc_ids;
  return;
}

void GRINS::DBCContainer::set_func( std::tr1::shared_ptr<libMesh::FunctionBase<Number> > func )
{
  _func = func;
  return;
}

std::vector<GRINS::VariableName> GRINS::DBCContainer::get_var_names() const
{
  return _var_names;
}

std::set<GRINS::BoundaryID> GRINS::DBCContainer::get_bc_ids() const
{
  return _bc_ids;
}

std::tr1::shared_ptr<libMesh::FunctionBase<Number> > GRINS::DBCContainer::get_func() const
{
  return _func;
}
