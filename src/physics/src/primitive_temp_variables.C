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

// This class
#include "grins/primitive_temp_variables.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

// GRINS
#include "grins/variable_name_defaults.h"

namespace GRINS
{
  PrimitiveTempVariables::PrimitiveTempVariables( const GetPot& input )
    : _T_var_name( input("Physics/VariableNames/Temperature", T_var_name_default ) )
  {
    return;
  }

  PrimitiveTempVariables::~PrimitiveTempVariables()
  {
    return;
  }

  void PrimitiveTempVariables::init( libMesh::FEMSystem* system )
  {
    libmesh_assert( system->has_variable(_T_var_name) );

    _T_var = system->variable_number( _T_var_name );
    
    return;
  }

} // end namespace GRINS
