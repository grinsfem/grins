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


#ifndef GRINS_STEADY_SOLVER_H
#define GRINS_STEADY_SOLVER_H

//GRINS
#include "grins/solver.h"

namespace GRINS
{
  class SteadySolver : public Solver
  {
  public:

    using Solver::Solver;

    virtual ~SteadySolver() = default;

    virtual void solve( SolverContext& context ) override;

    virtual void adjoint_qoi_parameter_sensitivity
    (SolverContext&                  context,
     const libMesh::QoISet&          qoi_indices,
     const libMesh::ParameterVector& parameters_in,
     libMesh::SensitivityData&       sensitivities) const override;

    virtual void forward_qoi_parameter_sensitivity
    (SolverContext&                  context,
     const libMesh::QoISet&          qoi_indices,
     const libMesh::ParameterVector& parameters_in,
     libMesh::SensitivityData&       sensitivities) const override;

  protected:

    virtual void init_time_solver(MultiphysicsSystem* system) override;

  };
} // namespace GRINS
#endif // GRINS_STEADY_SOLVER_H
