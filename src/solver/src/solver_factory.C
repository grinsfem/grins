//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
#include "grins/solver_factory.h"

// GRINS
#include "grins/grins_steady_solver.h"
#include "grins/grins_unsteady_solver.h"
#include "grins/steady_mesh_adaptive_solver.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  SolverFactory::SolverFactory()
  {
    return;
  }

  SolverFactory::~SolverFactory()
  {
    return;
  }

  std::tr1::shared_ptr<Solver> SolverFactory::build(const GetPot& input)
  {
    bool mesh_adaptive = input("MeshAdaptivity/mesh_adaptive", false );

    bool transient = input("unsteady-solver/transient", false );

    std::tr1::shared_ptr<Solver> solver;  // Effectively NULL

    if(transient && !mesh_adaptive)
      {
        solver.reset( new UnsteadySolver(input) );
      }
    else if( !transient && !mesh_adaptive )
      {
        solver.reset( new SteadySolver(input) );
      }
    else if( !transient && mesh_adaptive )
      {
        solver.reset( new SteadyMeshAdaptiveSolver(input) );
      }
    else if( transient && mesh_adaptive )
      {
        libmesh_not_implemented();
      }
    else
      {
        std::cerr << "Invalid solver options!" << std::endl;
      }

    return solver;
  }

} // namespace GRINS
