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
#include "grins/primitive_flow_variables.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

// GRINS
#include "grins/variable_name_defaults.h"

namespace GRINS
{
  PrimitiveFlowVariables::PrimitiveFlowVariables( const GetPot& input )
    :  _u_var_name( input("Physics/VariableNames/u_velocity", u_var_name_default ) ),
       _v_var_name( input("Physics/VariableNames/v_velocity", v_var_name_default ) ),
       _w_var_name( input("Physics/VariableNames/w_velocity", w_var_name_default ) ),
       _p_var_name( input("Physics/VariableNames/pressure",   p_var_name_default ) )
  {
    return;
  }

  PrimitiveFlowVariables::~PrimitiveFlowVariables()
  {
    return;
  }

  void PrimitiveFlowVariables::init( libMesh::FEMSystem* system )
  {
    libmesh_assert( system->has_variable( _u_var_name ) );
    libmesh_assert( system->has_variable( _v_var_name ) );
    libmesh_assert( system->has_variable( _p_var_name ) );

    _u_var = system->variable_number( _u_var_name );
    _v_var = system->variable_number( _v_var_name );

    if ( system->get_mesh().mesh_dimension() == 3)
      {
        libmesh_assert( system->has_variable( _w_var_name ) );
        _w_var = system->variable_number( _w_var_name );
      }

    _p_var = system->variable_number( _p_var_name );

    return;
  }

} // end namespace GRINS
