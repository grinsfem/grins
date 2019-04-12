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

#ifndef GRINS_NEUMANN_BC_CONTAINER_H
#define GRINS_NEUMANN_BC_CONTAINER_H

// GRINS
#include "grins/var_typedefs.h"
#include "grins/neumann_bc_abstract.h"
#include <memory>
#include "grins/fe_variables_base.h"

// libMesh
#include "libmesh/auto_ptr.h" // UniquePtr

namespace GRINS
{
  class NeumannBCContainer
  {
  public:
    NeumannBCContainer( const std::set<BoundaryID>& bc_ids, const FEVariablesBase& fe_var,
                        std::shared_ptr<NeumannBCAbstract>& func )
      : _bc_ids(bc_ids),
        _fe_var(fe_var),
        _func(func)
    {}

    ~NeumannBCContainer(){};

    bool has_bc_id( BoundaryID bc_id )
    { return (_bc_ids.find(bc_id) != _bc_ids.end()); }

    const FEVariablesBase& get_fe_var()
    { return _fe_var; }

    std::shared_ptr<NeumannBCAbstract>& get_func()
    { return _func; }

  private:

    std::set<BoundaryID> _bc_ids;

    const FEVariablesBase& _fe_var;

    std::shared_ptr<NeumannBCAbstract> _func;
  };

} // end namespace GRINS

#endif // GRINS_NEUMANN_BC_CONTAINER_H
