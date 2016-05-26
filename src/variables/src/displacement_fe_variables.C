//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  DisplacementFEVariables::DisplacementFEVariables( const GetPot& input,
                                                    const std::string& physics_name,
                                                    bool is_2D, bool is_3D,
                                                    bool is_constraint_var )
    :  MultiVarSingleFETypeVariable(input,physics_name,"",this->old_var_names(),this->default_names(),
                                    this->subsection(),"LAGRANGE","FIRST",is_constraint_var),
       _have_v(false),
       _have_w(false),
       _u_idx(0),
       _v_idx(1),
       _w_idx(2),
       _is_2D(is_2D),
       _is_3D(is_3D)
  {}

  void DisplacementFEVariables::init( libMesh::FEMSystem* system )
  {
    unsigned int mesh_dim = system->get_mesh().mesh_dimension();

    // The order matters here. We *must* do w first since we use pop_back().
    if ( mesh_dim == 3 || _is_3D )
      _have_w = true;
    else
      _var_names.pop_back();

    if ( mesh_dim >= 2 || _is_2D || _is_3D)
        _have_v = true;
    else
        _var_names.pop_back();

    SingleFETypeVariable::init(system);
  }

} // end namespace GRINS
