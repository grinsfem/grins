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


#ifndef GRINS_UNSTEADY_SOLVER_H
#define GRINS_UNSTEADY_SOLVER_H

//GRINS
#include "grins/grins_solver.h"

//libMesh
#include "libmesh/system_norm.h"

namespace GRINS
{
  class UnsteadySolver : public Solver
  {
  public:

    UnsteadySolver( const GetPot& input );
    virtual ~UnsteadySolver();

    virtual void solve( SolverContext& context );

  protected:

    virtual void init_time_solver(GRINS::MultiphysicsSystem* system);

    unsigned int _n_timesteps;
    unsigned int _backtrack_deltat;
    double _theta;
    double _deltat;

    // Options for adaptive time solvers
    double _target_tolerance;
    double _upper_tolerance;
    double _max_growth;
    SystemNorm _component_norm;
  };

} // end namespace GRINS
#endif // GRINS_UNSTEADY_SOLVER_H
