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

#ifndef GRINS_DISPLACEMENT_CONTINUATION_SOLVER_FACTORY_H
#define GRINS_DISPLACEMENT_CONTINUATION_SOLVER_FACTORY_H

// GRINS
#include "grins/solver_factory.h"
#include "grins/solver_parsing.h"

#include "displacement_continuation_solver.h"

namespace GRINS
{
  class DisplacementContinuationSolverFactory : public SolverFactory
  {
  public:
    DisplacementContinuationSolverFactory(){};
    virtual ~DisplacementContinuationSolverFactory(){};

    virtual SharedPtr<GRINS::Solver> build(const GetPot& input);
  };

  SharedPtr<GRINS::Solver> DisplacementContinuationSolverFactory::build(const GetPot& input)
  {
     std::cout << "HELLO!!!" << std::endl;

    std::string solver_type = input("SolverOptions/solver_type", "DIE!");

    SharedPtr<Solver> solver;  // Effectively NULL

    if( solver_type == std::string("displacement_continuation") )
      {
        solver.reset( new DisplacementContinuationSolver(input) );
      }
    else
      solver = SolverFactory::build(input);

    return solver;
  }

} // end namespace GRINS

#endif // GRINS_DISPLACEMENT_CONTINUATION_SOLVER_FACTORY_H
