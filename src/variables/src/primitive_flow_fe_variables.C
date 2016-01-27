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
#include "grins/primitive_flow_fe_variables.h"

// GRINS
#include "grins/variable_name_defaults.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  PrimitiveFlowFEVariables::PrimitiveFlowFEVariables( const GetPot& input, const std::string& physics_name )
    :  FEVariablesBase(),
       PrimitiveFlowVariables(input),
       _u_fe_idx(0),
       _p_fe_idx(1)
  {
    _family.resize(2,libMesh::INVALID_FE);
    _order.resize(2,libMesh::INVALID_ORDER);

    _family[_u_fe_idx] = libMesh::Utility::string_to_enum<GRINSEnums::FEFamily>( input("Physics/"+physics_name+"/V_FE_family", input("Physics/"+physics_name+"/FE_family", "LAGRANGE") ) );

    _family[_p_fe_idx ] = libMesh::Utility::string_to_enum<GRINSEnums::FEFamily>( input("Physics/"+physics_name+"/P_FE_family", input("Physics/"+physics_name+"/FE_family", "LAGRANGE") ) );

    _order[_u_fe_idx ] = libMesh::Utility::string_to_enum<GRINSEnums::Order>( input("Physics/"+physics_name+"/V_order", "SECOND") );

    _order[_p_fe_idx] = libMesh::Utility::string_to_enum<GRINSEnums::Order>( input("Physics/"+physics_name+"/P_order", "FIRST") );
  }

  void PrimitiveFlowFEVariables::init( libMesh::FEMSystem* system )
  {
    _vars[_u_idx] = system->add_variable( _var_names[_u_idx], _order[_u_fe_idx], _family[_u_fe_idx]);
    _vars[_v_idx] = system->add_variable( _var_names[_v_idx], _order[_u_fe_idx], _family[_u_fe_idx]);

    if ( system->get_mesh().mesh_dimension() == 3)
      _vars[_w_idx] = system->add_variable( _var_names[_w_idx], _order[_u_fe_idx], _family[_u_fe_idx]);

    _vars[_p_idx] = system->add_variable( _var_names[_p_idx], _order[_p_fe_idx], _family[_p_fe_idx]);
  }

} // end namespace GRINS
