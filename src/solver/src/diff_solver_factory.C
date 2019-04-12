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

#include "grins/diff_solver_factory.h"

namespace GRINS
{
  std::unique_ptr<libMesh::DiffSolver> DiffSolverFactoryAbstract::create()
  {
    if(!_system)
      libmesh_error_msg("ERROR: Must must MultiphysicsSystem pointer before calling DiffSolverFactoryAbstract::create()!");

    return this->build_diff_solver(*_system);
  }

  // Full specialization for the Factory<libMesh::DiffSolver>
  template<>
  std::map<std::string, FactoryAbstract<libMesh::DiffSolver>*>&
  FactoryAbstract<libMesh::DiffSolver>::factory_map()
  {
    static std::map<std::string, FactoryAbstract<libMesh::DiffSolver>*> _map;
    return _map;
  }

  // Definition of static members
  MultiphysicsSystem * DiffSolverFactoryAbstract::_system = NULL;

} // end namespace GRINS
