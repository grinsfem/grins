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
#include "grins/primitive_flow_variables.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

// GRINS
#include "grins/variable_name_defaults.h"

namespace GRINS
{
  PrimitiveFlowVariables::PrimitiveFlowVariables( const GetPot& input )
    :  VariablesBase(),
       _u_idx(0),
       _v_idx(1),
       _w_idx(2),
       _p_idx(3)
  {
    _vars.resize(4,invalid_var_index);
    _var_names.resize(4);

    _var_names[_u_idx] = input("Physics/VariableNames/u_velocity", u_var_name_default );
    _var_names[_v_idx] = input("Physics/VariableNames/v_velocity", v_var_name_default );
    _var_names[_w_idx] = input("Physics/VariableNames/w_velocity", w_var_name_default );
    _var_names[_p_idx] = input("Physics/VariableNames/pressure", p_var_name_default );
  }

  void PrimitiveFlowVariables::init( libMesh::FEMSystem* system )
  {
    libmesh_assert( system->has_variable( _var_names[_u_idx] ) );
    libmesh_assert( system->has_variable( _var_names[_v_idx] ) );
    libmesh_assert( system->has_variable( _var_names[_p_idx] ) );

    _vars[_u_idx] = system->variable_number( _var_names[_u_idx] );
    _vars[_v_idx] = system->variable_number( _var_names[_v_idx] );

    if ( system->get_mesh().mesh_dimension() == 3)
      {
        libmesh_assert( system->has_variable( _var_names[_w_idx] ) );
        _vars[_w_idx] = system->variable_number( _var_names[_w_idx] );
      }

    _vars[_p_idx] = system->variable_number( _var_names[_p_idx] );
  }

} // end namespace GRINS
