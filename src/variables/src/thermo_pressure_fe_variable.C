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
#include "grins/thermo_pressure_fe_variable.h"

// GRINS
#include "grins/grins_enums.h"
#include "grins/variable_name_defaults.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  ThermoPressureFEVariable::ThermoPressureFEVariable( const GetPot& input, const std::string& /*physics_name*/ )
    :  ThermoPressureVariable(input),
       _P_FE_family( libMesh::SCALAR ),
       _P_order( libMesh::FIRST )
  {}

  void ThermoPressureFEVariable::init( libMesh::FEMSystem* system )
  {
    _vars[0] = system->add_variable( _var_names[0], libMesh::FIRST, libMesh::SCALAR );
  }

} // end namespace GRINS
