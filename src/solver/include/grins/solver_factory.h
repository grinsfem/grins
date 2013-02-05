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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_SOLVER_FACTORY_H
#define GRINS_SOLVER_FACTORY_H

//GRINS
#include "grins/grins_solver.h"

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  //! This object handles constructing the solver to be used.
  /*! To allow the user to easily extend the (limited) available solvers,
      the solver construction is handled in this object. */
  class SolverFactory
  {
  public:

    SolverFactory();
    virtual ~SolverFactory();

    //! Builds GRINS::Solver object.
    /*! Users should override this method to construct 
        their own solvers. */
    virtual std::tr1::shared_ptr<GRINS::Solver> build(const GetPot& input);

  };
} // namespace GRINS
#endif //GRINS_SOLVER_FACTORY_H
