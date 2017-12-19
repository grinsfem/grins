//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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

#include "grins/shared_ptr.h"

// GRINS
#include "grins/shared_ptr.h"
#include "grins/visualization.h"
#include "grins/postprocessed_quantities.h"
#include "grins/qoi_output.h"

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
    ~SolverContext(){};

    GRINS::MultiphysicsSystem* system;
    std::shared_ptr<libMesh::EquationSystems> equation_system;
    std::shared_ptr<GRINS::Visualization> vis;
    unsigned int timesteps_per_vis;
    unsigned int timesteps_per_perflog;
    bool output_vis;
    bool output_adjoint;
    bool output_residual;
    bool output_residual_sensitivities;
    bool output_solution_sensitivities;
    bool print_perflog;
    bool print_scalars;
    bool do_adjoint_solve;

    std::shared_ptr<QoIOutput> qoi_output;

    std::shared_ptr<PostProcessedQuantities<libMesh::Real> > postprocessing;

    std::shared_ptr<libMesh::ErrorEstimator> error_estimator;

    bool have_restart;

  };

} // end namespace GRINS
#endif // GRINS_SOLVER_CONTEXT_H
