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

#include "grins/nbc_container.h"

namespace GRINS
{
  NBCContainer::NBCContainer()
  {
    return;
  }

  NBCContainer::~NBCContainer()
  {
    return;
  }

  void NBCContainer::set_bc_id( BoundaryID bc_id )
  {
    _bc_id = bc_id;
    return;
  }

  BoundaryID NBCContainer::get_bc_id() const
  {
    return _bc_id;
  }

  void NBCContainer::add_var_func_pair( VariableIndex var, 
					std::tr1::shared_ptr<NeumannFuncObj> func )
  {
    if( _funcs.find(var) != _funcs.end() )
      {
	std::cerr << "Error: Can only specify one function per variable" << std::endl;
	libmesh_error();
      }

    _funcs.insert( std::make_pair( var, func ) );
    return;
  }

  std::tr1::shared_ptr<NeumannFuncObj> NBCContainer::get_func( VariableIndex var ) const
  {
    libmesh_assert( _funcs.find(var) != _funcs.end() );
    return _funcs.find(var)->second;
  }

} // namespace GRINS
