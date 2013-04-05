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
#include "grins/neumann_func_obj.h"

// GRINS
#include "grins/cached_values.h"

// libMesh
#include "libmesh/fem_context.h"

namespace GRINS
{

  NeumannFuncObj::NeumannFuncObj( )
  {
    return;
  }

  NeumannFuncObj::~NeumannFuncObj( )
  {
    return;
  }

  libMesh::Point NeumannFuncObj::value( const libMesh::FEMContext&,
					const CachedValues&,
					const unsigned int )
  {
    // By default, does nothing.
    /* \todo Should we libmesh_error() instead?*/
    return libMesh::Point();
  }

  libMesh::Real NeumannFuncObj::normal_value( const libMesh::FEMContext&,
					      const CachedValues&,
					      const unsigned int )
  {
    // By default, does nothing.
    /* \todo Should we libmesh_error() instead?*/
    return 0.0;
  }

  libMesh::Point NeumannFuncObj::derivative( const libMesh::FEMContext&,
					     const CachedValues&,
					     const unsigned int )
  {
    // By default, does nothing.
    /* \todo Should we libmesh_error() instead?*/
    return libMesh::Point(0.0,0.0,0.0);
  }

  libMesh::Point NeumannFuncObj::derivative( const libMesh::FEMContext&,
					     const CachedValues&,
					     const unsigned int,
					     const VariableIndex )
  {
    // By default, does nothing.
    /* \todo Should we libmesh_error() instead?*/
    return libMesh::Point(0.0,0.0,0.0);
  }
  

  libMesh::Real NeumannFuncObj::normal_derivative( const libMesh::FEMContext&,
						   const CachedValues&,
						   const unsigned )
  {
    // By default, does nothing.
    /* \todo Should we libmesh_error() instead?*/
    return 0.0;
  }

  libMesh::Real NeumannFuncObj::normal_derivative( const libMesh::FEMContext&,
						   const CachedValues&,
						   const unsigned int, 
						   const VariableIndex )
  {
    // By default, does nothing.
    /* \todo Should we libmesh_error() instead?*/
    return 0.0;
  }

  const std::vector<VariableIndex>& NeumannFuncObj::get_other_jac_vars()
  {
    return _jac_vars;
  }

} // namespace GRINS

