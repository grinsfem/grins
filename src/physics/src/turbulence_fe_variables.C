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
#include "grins/turbulence_fe_variables.h"

// GRINS
#include "grins/grins_enums.h"
#include "grins/variable_name_defaults.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  TurbulenceFEVariables::TurbulenceFEVariables( const GetPot& input, const std::string& physics_name )
    :  TurbulenceVariables(input),
       _TU_FE_family( libMesh::Utility::string_to_enum<GRINSEnums::FEFamily>( input("Physics/"+physics_name+"/TU_FE_family", input("Physics/"+physics_name+"/FE_family", "LAGRANGE") ) ) ),
       _TU_order( libMesh::Utility::string_to_enum<GRINSEnums::Order>( input("Physics/"+physics_name+"/TU_order", "FIRST") ) )
  {
    return;
  }

  TurbulenceFEVariables::~TurbulenceFEVariables()
  {
    return;
  }

  void TurbulenceFEVariables::init( libMesh::FEMSystem* system )
  {
    _nu_var = system->add_variable( _nu_var_name, this->_TU_order, _TU_FE_family);     
    return;
  }

} // end namespace GRINS
