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


#ifndef GRINS_DBC_CONTAINER_H
#define GRINS_DBC_CONTAINER_H

#include <vector>
#include <set>

#include "boost/tr1/memory.hpp"

// libMesh
#include "libmesh/fem_function_base.h"
#include "libmesh/function_base.h"

// GRINS
#include "grins/var_typedefs.h"

namespace GRINS
{
  //! Simple helper class to setup general Dirichlet boundary conditions
  /*! This class is to temporarily stash data necessary for setting
      up libMesh::DirichletBoundary objects. Actual instantiation
      of libMesh::DirichletBoundary object is handled internally by
      GRINS::Physics::init_user_dirichlet_bcs. For each Dirichlet bc
      function, there is a unique DBCContainer object. */
  class DBCContainer
  {
  public:
    
    DBCContainer();
    ~DBCContainer();

    //! Add variables that are constrained by the Dirichlet bc
    void add_var_name( const GRINS::VariableName& var );

    //! Add boundary id's for which this Dirichlet bc is to be applied
    void add_bc_id( const GRINS::BoundaryID bc_id );

    //! Add the Dirichlet bc functor
    /*! There is only one Dirichlet bc functor for each container object */
    void set_func( std::tr1::shared_ptr<libMesh::FunctionBase<libMesh::Number> > func );

    void set_fem_func_string( const std::string& s );

    std::vector<GRINS::VariableName> get_var_names() const;
    std::set<GRINS::BoundaryID> get_bc_ids() const;

    std::tr1::shared_ptr<libMesh::FunctionBase<libMesh::Number> >
    get_func() const;

    const std::string&
    get_fem_func_string() const;

  private:

    std::vector<GRINS::VariableName> _var_names;
    std::set<GRINS::BoundaryID> _bc_ids;
    std::tr1::shared_ptr<libMesh::FunctionBase<libMesh::Number> > _func;
    std::string _parsed_fem_func_string;
    
  };
}
#endif // GRINS_DBC_CONTAINER_H
