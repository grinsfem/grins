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
#include "grins/displacement_fe_variables.h"

// GRINS
#include "grins/variable_name_defaults.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  DisplacementFEVariables::DisplacementFEVariables( const GetPot& input, const std::string& physics_name )
    :  FEVariablesBase(),
       DisplacementVariables(input)
  {
    _family.resize(1, libMesh::Utility::string_to_enum<GRINSEnums::FEFamily>( input("Physics/"+physics_name+"/FE_family", "LAGRANGE") ) );
    _order.resize(1, libMesh::Utility::string_to_enum<GRINSEnums::Order>( input("Physics/"+physics_name+"/order", "FIRST") ) );
  }

  void DisplacementFEVariables::init( libMesh::FEMSystem* system, bool is_2D, bool is_3D )
  {
    _vars[_u_idx] = system->add_variable( _var_names[_u_idx], this->_order[0], _family[0]);

    if ( system->get_mesh().mesh_dimension() >= 2 || is_2D || is_3D)
      {
        _have_v = true;
        _vars[_v_idx] = system->add_variable( _var_names[_v_idx], this->_order[0], _family[0]);
      }

    if ( system->get_mesh().mesh_dimension() == 3 || is_3D )
      {
        _have_w = true;
        _vars[_w_idx] = system->add_variable( _var_names[_w_idx], this->_order[0], _family[0]);
      }
  }

} // end namespace GRINS
