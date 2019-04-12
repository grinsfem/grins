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

#ifndef GRINS_SOLVER_NAMES_H
#define GRINS_SOLVER_NAMES_H

// C++
#include <string>

namespace GRINS
{
  class SolverNames
  {
  public:

    static const std::string steady_solver()
    { return "grins_steady_solver"; }

    static const std::string unsteady_solver()
    { return "grins_unsteady_solver"; }

    static const std::string steady_mesh_adaptive_solver()
    { return "grins_steady_mesh_adaptive_solver"; }

    static const std::string unsteady_mesh_adaptive_solver()
    { return "grins_unsteady_mesh_adaptive_solver"; }

    static const std::string libmesh_euler_solver()
    { return "libmesh_euler_solver"; }

    static const std::string libmesh_euler2_solver()
    { return "libmesh_euler2_solver"; }

    static const std::string libmesh_newmark_solver()
    { return "libmesh_newmark"; }

  };
} // end namespace GRINS

#endif // GRINS_SOLVER_NAMES_H
