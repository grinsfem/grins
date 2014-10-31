//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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
#include "grins/solid_mechanics_variables.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

// GRINS
#include "grins/variable_name_defaults.h"

namespace GRINS
{
  SolidMechanicsVariables::SolidMechanicsVariables( const GetPot& input )
    : _u_var_name( input("Physics/VariableNames/u_displacment", u_disp_name_default ) ),
      _v_var_name( input("Physics/VariableNames/v_displacment", v_disp_name_default ) ),
      _w_var_name( input("Physics/VariableNames/w_displacment", w_disp_name_default ) )
  {
    return;
  }

  SolidMechanicsVariables::~SolidMechanicsVariables()
  {
    return;
  }

  void SolidMechanicsVariables::init( libMesh::FEMSystem* system, bool is_2D, bool is_3D )
  {
    libmesh_assert( system->has_variable( _u_var_name ) );
    _u_var = system->variable_number( _u_var_name );

    if ( system->get_mesh().mesh_dimension() >= 2 || is_2D )
      {
        libmesh_assert( system->has_variable( _v_var_name ) );
        _v_var = system->variable_number( _v_var_name );
      }

    if ( system->get_mesh().mesh_dimension() == 3 || is_3D )
      {
        libmesh_assert( system->has_variable( _w_var_name ) );
        _w_var = system->variable_number( _w_var_name );
      }

    return;
  }

} // end namespace GRINS
