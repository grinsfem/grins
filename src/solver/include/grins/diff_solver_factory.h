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

#ifndef GRINS_DIFF_SOLVER_FACTORY_H
#define GRINS_DIFF_SOLVER_FACTORY_H

// GRINS
#include "grins/factory_abstract.h"
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/diff_solver.h"

namespace GRINS
{
  class DiffSolverFactoryAbstract : public FactoryAbstract<libMesh::DiffSolver>
  {
  public:
    DiffSolverFactoryAbstract( const std::string & diff_solver_name )
      : FactoryAbstract<libMesh::DiffSolver>(diff_solver_name)
    {}

    virtual ~DiffSolverFactoryAbstract() =0;

    static void set_system( MultiphysicsSystem * system )
    { _system = system; }

  protected:

   static MultiphysicsSystem * _system;

   virtual std::unique_ptr<libMesh::DiffSolver> build_diff_solver( MultiphysicsSystem & system ) =0;

  private:

    virtual std::unique_ptr<libMesh::DiffSolver> create();

    DiffSolverFactoryAbstract();
  };

  inline
  DiffSolverFactoryAbstract::~DiffSolverFactoryAbstract(){}


  template<typename DiffSolverType>
  class DiffSolverFactoryBasic : public DiffSolverFactoryAbstract
  {
  public:

    DiffSolverFactoryBasic( const std::string & diff_solver_name )
      : DiffSolverFactoryAbstract(diff_solver_name)
    {}

    virtual ~DiffSolverFactoryBasic(){}

  protected:

    virtual std::unique_ptr<libMesh::DiffSolver> build_diff_solver( MultiphysicsSystem & system )
    { return std::unique_ptr<libMesh::DiffSolver>( new DiffSolverType(system) ); }

  private:

    DiffSolverFactoryBasic();

  };

} // end namespace GRINS

#endif // GRINS_DIFF_SOLVER_FACTORY_H
