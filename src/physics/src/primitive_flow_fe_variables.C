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
#include "grins/primitive_flow_fe_variables.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_system.h"

// GRINS
#include "grins/variable_name_defaults.h"

namespace GRINS
{
  PrimitiveFlowFEVariables::PrimitiveFlowFEVariables( const GetPot& input, const std::string& physics_name )
    :  PrimitiveFlowVariables(input),
       _V_FE_family( libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+physics_name+"/V_FE_family", input("Physics/"+physics_name+"/FE_family", "LAGRANGE") ) ) ),
       _P_FE_family( libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+physics_name+"/P_FE_family", input("Physics/"+physics_name+"/FE_family", "LAGRANGE") ) ) ),
       _V_order( libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+physics_name+"/V_order", "SECOND") ) ),
       _P_order( libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+physics_name+"/P_order", "FIRST") ) )
  {
    return;
  }

  PrimitiveFlowFEVariables::~PrimitiveFlowFEVariables()
  {
    return;
  }

  void PrimitiveFlowFEVariables::init( libMesh::FEMSystem* system )
  {
    _u_var = system->add_variable( _u_var_name, this->_V_order, _V_FE_family);
    _v_var = system->add_variable( _v_var_name, this->_V_order, _V_FE_family);

    if ( system->get_mesh().mesh_dimension() == 3)
      _w_var = system->add_variable( _w_var_name, this->_V_order, _V_FE_family);

    _p_var = system->add_variable( _p_var_name, this->_P_order, _P_FE_family);

    return;
  }

} // end namespace GRINS
