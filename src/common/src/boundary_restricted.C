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
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/getpot.h"

// C++
#include <algorithm>

namespace GRINS
{
  void BoundaryRestricted::init_bcids
  (const GetPot& input,
   std::string_view name)
  {
    // Are we computing on boundary sides identified by nodes? (false by default)
    const std::string nodal_sides_var = std::string(name)+"/nodal_sides";
    _nodal_sides = input(nodal_sides_var, false);

    // Read boundary ids for which we want to compute
    const std::string bc_ids_var = std::string(name)+"/bc_ids";
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
    if (!context.has_elem())
      return false;

    for (auto id : _bc_ids)
      if( context.has_side_boundary_id(id) )
        return true;

    if (_nodal_sides)
      {
        const libMesh::BoundaryInfo & boundary_info =
          context.get_multiphysics_system().get_mesh().get_boundary_info();
        const libMesh::Elem & elem = context.get_elem();
        const unsigned char s = context.get_side();
        auto side_nodes = elem.nodes_on_side(s);
        for (auto id : _bc_ids)
          {
            bool all_have_id = true;
            for (auto n : side_nodes)
              if (!boundary_info.has_boundary_id(elem.node_ptr(n), id))
                {
                  all_have_id = false;
                  break;
                }
            if (all_have_id)
              return true;
          }
      }

    return false;
  }
} //namespace GRINS
