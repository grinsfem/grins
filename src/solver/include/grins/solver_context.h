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


#ifndef GRINS_SOLVER_CONTEXT_H
#define GRINS_SOLVER_CONTEXT_H

#include "boost/tr1/memory.hpp"

// GRINS
#include "grins/visualization.h"
#include "grins/postprocessed_quantities.h"

// libMesh
#include "libmesh/error_estimator.h"
#include "libmesh/equation_systems.h"

namespace GRINS
{
  // Forward declarations
  class MultiphysicsSystem;

  //! Simple class to hold objects passed to Solver::solve
  /*! Allows some flexibility for adding objects needed to pass to the Solver::solve
      method so that the solver can still be agnostic to creation etc. of those objects,
      but can operate on them. 
   */
  class SolverContext
  {
  public:
    
    SolverContext();
    ~SolverContext();

    GRINS::MultiphysicsSystem* system;
    std::tr1::shared_ptr<libMesh::EquationSystems> equation_system;
    std::tr1::shared_ptr<GRINS::Visualization> vis;
    unsigned int timesteps_per_vis;
    bool output_vis;
    bool output_residual;

    std::tr1::shared_ptr<PostProcessedQuantities<Real> > postprocessing;

    std::tr1::shared_ptr<libMesh::ErrorEstimator> error_estimator;

  };

} // end namespace GRINS
#endif // GRINS_SOLVER_CONTEXT_H
