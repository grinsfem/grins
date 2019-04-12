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
#include "grins/solver_factory_initializer.h"

// GRINS
#include "grins/solver_factory_basic.h"
#include "grins/solver_names.h"

// GRINS-Solvers
#include "grins/steady_solver.h"
#include "grins/unsteady_solver.h"
#include "grins/steady_mesh_adaptive_solver.h"
#include "grins/unsteady_mesh_adaptive_solver.h"

namespace GRINS
{
  SolverFactoryInitializer::SolverFactoryInitializer()
  {
    static SolverFactoryBasic<UnsteadySolver>
      grins_factory_unsteady_solver(SolverNames::unsteady_solver());

    static SolverFactoryBasic<SteadySolver>
      grins_factory_steady_solver(SolverNames::steady_solver());

    static SolverFactoryBasic<UnsteadyMeshAdaptiveSolver>
      grins_factory_unsteady_mesh_adapt_solver(SolverNames::unsteady_mesh_adaptive_solver());

    static SolverFactoryBasic<SteadyMeshAdaptiveSolver>
      grins_factory_steady_mesh_adapt_solver(SolverNames::steady_mesh_adaptive_solver());
  }
} // end namespace GRINS
