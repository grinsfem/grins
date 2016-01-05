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
#include "grins/thermo_pressure_variable.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  ThermoPressureVariable::ThermoPressureVariable( const GetPot& input )
    : VariablesBase()
  {
    _vars.resize(1,invalid_var_index);
    _var_names.resize(1, input("Physics/VariableNames/thermo_presure", "p0" ) );
  }

  void ThermoPressureVariable::init( libMesh::FEMSystem* system )
  {
    libmesh_assert( system->has_variable( _var_names[0] ) );

    _vars[0] = system->variable_number( _var_names[0] );
  }

} // end namespace GRINS
