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
#include "grins/displacement_variables.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

// GRINS
#include "grins/variable_name_defaults.h"

namespace GRINS
{
  DisplacementVariables::DisplacementVariables( const GetPot& input )
    : VariablesBase(),
      _have_v(false),
      _have_w(false),
       _u_idx(0),
       _v_idx(1),
       _w_idx(2)
  {
    _vars.resize(3,invalid_var_index);
    _var_names.resize(3);

    _var_names[_u_idx] = input("Physics/VariableNames/u_displacment", u_disp_name_default );
    _var_names[_v_idx] = input("Physics/VariableNames/v_displacment", v_disp_name_default );
    _var_names[_w_idx] = input("Physics/VariableNames/w_displacment", w_disp_name_default );
  }

  void DisplacementVariables::init( libMesh::FEMSystem* system )
  {
    libmesh_assert( system->has_variable( _var_names[_u_idx] ) );
    _vars[_u_idx] = system->variable_number( _var_names[_u_idx] );

    if ( system->has_variable( _var_names[_v_idx] ) )
      {
        _have_v = true;
        _vars[_v_idx] = system->variable_number( _var_names[_v_idx] );
      }

    if ( system->has_variable( _var_names[_w_idx] ) )
      {
        _have_w = true;
        _vars[_w_idx] = system->variable_number( _var_names[_w_idx] );
      }

    return;
  }

} // end namespace GRINS
