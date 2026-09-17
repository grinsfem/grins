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


// This class
#include "grins/boundary_restricted.h"

// GRINS
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/getpot.h"

// C++
#include <algorithm>

namespace GRINS
{
  void BoundaryRestricted::init_bcids
  (const GetPot& input,
   std::string_view name)
  {
    // Read boundary ids for which we want to compute
    const std::string bc_ids_var = "QoI/"+std::string(name)+"/bc_ids";
    int num_bcs =  input.vector_variable_size(bc_ids_var);

    if( num_bcs <= 0 )
      {
        std::cerr << "Error: Must specify at least one boundary id to compute "
                  << name << std::endl
                  << "Found: " << num_bcs << std::endl;
        libmesh_error();
      }

    for( int i = 0; i < num_bcs; i++ )
      _bc_ids.push_back( input(bc_ids_var, -1, i ) );

    // We use a vector for efficiency but we should still handle any
    // duplicates from the user
    std::sort( _bc_ids.begin(), _bc_ids.end() );
    _bc_ids.erase( std::unique( _bc_ids.begin(), _bc_ids.end() ), _bc_ids.end() );
  }

  bool BoundaryRestricted::is_on_active_boundary( const AssemblyContext& context )
  {
    for (auto id : _bc_ids)
      if( context.has_side_boundary_id(id) )
        return true;

    return false;
  }
} //namespace GRINS
