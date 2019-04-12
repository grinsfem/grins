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

#ifndef GRINS_SOLVER_FACTORY_ABSTRACT_H
#define GRINS_SOLVER_FACTORY_ABSTRACT_H

// GRINS
#include "grins/factory_with_getpot.h"
#include "grins/solver.h"

namespace GRINS
{
  // According to the standard, we need a declaration of the
  // specialization which precedes any automatic instantiation.
  template<> const GetPot* FactoryWithGetPot<Solver>::_input;

  class SolverFactoryAbstract : public FactoryWithGetPot<Solver>
  {
  public:
    SolverFactoryAbstract( const std::string & bc_type_name )
      : FactoryWithGetPot<Solver>(bc_type_name)
    {}

    virtual ~SolverFactoryAbstract() =0;

  protected:

    virtual std::unique_ptr<Solver> build_solver( const GetPot & input ) =0;

  private:

    virtual std::unique_ptr<Solver> create();

    SolverFactoryAbstract();
  };

  inline
  SolverFactoryAbstract:: ~SolverFactoryAbstract(){}

} // end namespace GRINS

#endif // GRINS_SOLVER_FACTORY_ABSTRACT_H
