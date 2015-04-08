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
#include "grins/turbulence_variables.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

// GRINS
#include "grins/variable_name_defaults.h"

namespace GRINS
{
  TurbulenceVariables::TurbulenceVariables( const GetPot& input )
    :  _nu_var_name( input("Physics/VariableNames/turbulent_viscosity", nu_var_name_default ) )
  {
    return;
  }

  TurbulenceVariables::~TurbulenceVariables()
  {
    return;
  }

  void TurbulenceVariables::init( libMesh::FEMSystem* system )
  {
    libmesh_assert( system->has_variable( _nu_var_name ) );
    
    _nu_var = system->variable_number( _nu_var_name );
    
    return;
  }

} // end namespace GRINS
