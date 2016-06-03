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
#include "grins/velocity_fe_variables.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  VelocityFEVariables::VelocityFEVariables( const GetPot& input, const std::string& physics_name,
                                            bool is_constraint_var)
    :  MultiVarSingleFETypeVariable(input,physics_name,"V_",this->old_var_names(),this->default_names(),
                                    this->subsection(),"LAGRANGE","SECOND",is_constraint_var),
       _u_idx(0),
       _v_idx(1),
       _w_idx(2)
  {}

  void VelocityFEVariables::init( libMesh::FEMSystem* system )
  {
    libmesh_assert_greater_equal(system->get_mesh().mesh_dimension(), 2);

    if ( system->get_mesh().mesh_dimension() < 3)
      _var_names.pop_back();

    SingleFETypeVariable::init(system);
  }

} // end namespace GRINS
