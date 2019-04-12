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
#include "grins/simulation_initializer.h"

// GRINS
#include "grins/error_estimator_factory_initializer.h"
#include "grins/physics_factory_initializer.h"
#include "grins/boundary_condition_factory_initializer.h"
#include "grins/variable_factory_initializer.h"
#include "grins/solver_factory_initializer.h"
#include "grins/diff_solver_factory_initializer.h"

namespace GRINS
{
  bool SimulationInitializer::_is_initialized = false;

  SimulationInitializer::SimulationInitializer()
  {
    if( !_is_initialized )
      {
        ErrorEstimatorFactoryInitializer error_est_init;
        PhysicsFactoryInitializer physics_init;
        BoundaryConditionFactoryInitializer bc_init;
        VariableFactoryInitializer var_init;
        SolverFactoryInitializer solver_init;
        DiffSolverFactoryInitializer diff_solver_init;

        _is_initialized = true;
      }
  }
} // end namespace GRINS
