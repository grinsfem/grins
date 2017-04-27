//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_SOLVER_FACTORY_BASIC_H
#define GRINS_SOLVER_FACTORY_BASIC_H

#include "grins/solver_factory_abstract.h"

namespace GRINS
{
  template<typename DerivedSolver>
  class SolverFactoryBasic : public SolverFactoryAbstract
  {
  public:
    SolverFactoryBasic( const std::string & bc_type_name )
      : SolverFactoryAbstract(bc_type_name)
    {}

    virtual ~SolverFactoryBasic(){}

  protected:

    virtual libMesh::UniquePtr<Solver> build_solver( const GetPot & input );

  private:

    SolverFactoryBasic();
  };

  template<typename DerivedSolver>
  inline
  libMesh::UniquePtr<Solver>
  SolverFactoryBasic<DerivedSolver>::build_solver( const GetPot & input )
  {
    return libMesh::UniquePtr<Solver>( new DerivedSolver(input) );
  }

} // end namespace GRINS

#endif // GRINS_SOLVER_FACTORY_BASIC_H
