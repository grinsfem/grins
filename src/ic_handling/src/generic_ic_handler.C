//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "grins/ic_handling_base.h"

// GRINS
#include "grins/generic_ic_handler.h"
#include "grins/string_utils.h"

// libMesh
#include "libmesh/fem_context.h"
#include "libmesh/fem_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/dof_map.h"
#include "libmesh/parsed_function.h"
#include "libmesh/const_function.h"

// C++
#include "sstream"

namespace GRINS
{
  GenericICHandler::GenericICHandler(const std::string& physics_name,
                                     const GetPot& input)
    : ICHandlingBase( physics_name )
  {
    std::string id_str = "Physics/"+_physics_name+"/ic_ids";
    std::string bc_str = "Physics/"+_physics_name+"/ic_types";
    std::string var_str = "Physics/"+_physics_name+"/ic_variables";
    std::string val_str = "Physics/"+_physics_name+"/ic_values";

    this->read_ic_data( input, id_str, bc_str, var_str, val_str );

    return;
  }

  GenericICHandler::~GenericICHandler()
  {
    return;
  }

} // namespace GRINS
