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

#ifndef GRINS_TESTING_UTILS_H
#define GRINS_TESTING_UTILS_H

// C++
#include <limits>

// libMesh
#include "libmesh/libmesh_common.h"

namespace GRINSTesting
{
  class TestingUtils
  {
  public:

    //! Get absolute tolerance from input relative tol
    /*! CppUnit uses absolute tolerance for double tests. This will
        let the user input a relative tolerance and an exact value
        to get back the corresponding absolute tolerance. */
    static libMesh::Real abs_tol_from_rel_tol( libMesh::Real exact, libMesh::Real rel_tol )
    {
      return std::abs(exact)*rel_tol;
    }

    //! Convenience function
    static libMesh::Real epsilon()
    {
      return std::numeric_limits<libMesh::Real>::epsilon();
    }
  };
}

#endif // GRINS_TESTING_UTILS_H
