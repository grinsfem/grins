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


#ifndef GRINS_BOUNDARY_RESTRICTED_H
#define GRINS_BOUNDARY_RESTRICTED_H

// GRINS forward declarations
namespace GRINS {
  class AssemblyContext;
}

// libMesh
#include "libmesh/id_types.h"

// libMesh forward declarations
class GetPot;

// C++
#include <vector>
#include <string_view>

namespace GRINS
{
  class BoundaryRestricted
  {
  public:
    virtual ~BoundaryRestricted() = default;

    void init_bcids( const GetPot& input,
                     std::string_view name );

    bool is_on_active_boundary( const AssemblyContext& context );

  private:
    //! List of boundaries on which we want to compute
    std::vector<libMesh::boundary_id_type> _bc_ids;
  };
}
#endif //GRINS_BOUNDARY_RESTRICTED_H
