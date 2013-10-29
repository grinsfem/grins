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
#include "grins/primitive_temp_fe_variables.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  PrimitiveTempFEVariables::PrimitiveTempFEVariables( const GetPot& input, const std::string& physics_name )
    : PrimitiveTempVariables(input),
      _T_FE_family( libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+physics_name+"/T_FE_family", input("Physics/"+physics_name+"/FE_family", "LAGRANGE") ) ) ),
      _T_order( libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+physics_name+"/T_order", "SECOND") ) )
  {
    return;
  }

  PrimitiveTempFEVariables::~PrimitiveTempFEVariables()
  {
    return;
  }

  void PrimitiveTempFEVariables::init( libMesh::FEMSystem* system )
  {
    _T_var = system->add_variable( _T_var_name, this->_T_order, _T_FE_family );
    
    return;
  }

} // end namespace GRINS
