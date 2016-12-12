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
#include "grins/constraint_builder.h"

// GRINS
#include "grins/constrained_points.h"
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/dof_map.h"

namespace GRINS
{
  libMesh::UniquePtr<libMesh::System::Constraint>
  ConstraintBuilder::build_constraint_object( const GetPot& input,
                                              MultiphysicsSystem& system )
  {
    // For now we only have the one constraint object option: we build
    // a ConstrainedPoints object, initialize it from the input
    // options, and set it to be the libMesh Constraint object.
    //
    // If we ever want to add more constraint object options, we'll
    // also need a MultiConstraint object to hold them.

    return libMesh::UniquePtr<libMesh::System::Constraint>
      (new ConstrainedPoints(input, system));
  }
} // end namespace GRINS
