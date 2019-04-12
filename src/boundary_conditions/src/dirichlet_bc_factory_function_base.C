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

// This class
#include "grins/dirichlet_bc_factory_function_base.h"

// libMesh
#include "libmesh/function_base.h"
#include "libmesh/fem_function_base.h"

namespace GRINS
{
  template<typename FunctionType>
  std::unique_ptr<libMesh::DirichletBoundary> DirichletBCFactoryFunctionBase<FunctionType>::create()
  {
    // Make sure all necessary state has been setup
    this->check_state();

    // We need to make a copy of var_names because subclasses may
    // modify what variable names are present to allow the boundary
    // condition factory implementation to apply only a subset of the
    // variables.
    std::vector<std::string> local_var_names = this->get_var_names();

    std::unique_ptr<FunctionType>
      func = this->build_func( *(this->_input), *(this->_system),
                               local_var_names, this->_section );

    libmesh_assert(func);

    std::vector<VariableIndex> local_var_indices;
    this->build_var_indices(*(this->_system), local_var_names, local_var_indices);

    std::unique_ptr<libMesh::DirichletBoundary> new_dbc =
      this->make_dirichlet_boundary( *(this->_bc_ids), *(this->_system),
                                     func, local_var_indices );

    // Reset state for error checking during next construction
    this->reset_state();

    return new_dbc;
  }

  template<>
  std::unique_ptr<libMesh::DirichletBoundary>
  DirichletBCFactoryFunctionBase<libMesh::FunctionBase<libMesh::Number> >::make_dirichlet_boundary
  ( const std::set<BoundaryID>& bc_ids,
    const libMesh::System& /*system*/,
    std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >& func,
    const std::vector<VariableIndex>& var_indices )
  {
    return std::unique_ptr<libMesh::DirichletBoundary>(new libMesh::DirichletBoundary(bc_ids, var_indices, func.get()));
  }

  template<>
  std::unique_ptr<libMesh::DirichletBoundary>
  DirichletBCFactoryFunctionBase<libMesh::FEMFunctionBase<libMesh::Number> >::make_dirichlet_boundary
  ( const std::set<BoundaryID>& bc_ids,
    const libMesh::System& system,
    std::unique_ptr<libMesh::FEMFunctionBase<libMesh::Number> >& func,
    const std::vector<VariableIndex>& var_indices )
  {
    return std::unique_ptr<libMesh::DirichletBoundary>( new libMesh::DirichletBoundary(bc_ids, var_indices, system, func.get()) );
  }
}

// Instantiate the factories
template class GRINS::DirichletBCFactoryFunctionBase<libMesh::FunctionBase<libMesh::Number> >;
template class GRINS::DirichletBCFactoryFunctionBase<libMesh::FEMFunctionBase<libMesh::Number> >;
