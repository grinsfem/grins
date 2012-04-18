//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef DBC_CONTAINER_H
#define DBC_CONTAINER_H

#include <vector>
#include <set>

#include "boost/tr1/memory.hpp"

// libMesh
#include "function_base.h"

// GRINS
#include "var_typedefs.h"

namespace GRINS
{
  class DBCContainer
  {
  public:
    
    DBCContainer();
    ~DBCContainer();

    void set_var_names( const std::vector<GRINS::VariableName>& var_names );
    void set_bc_ids( const std::set<GRINS::BoundaryID>& bc_ids );
    void set_func(std::tr1::shared_ptr<libMesh::FunctionBase<Number> > func );

    std::vector<GRINS::VariableName> get_var_names() const;
    std::set<GRINS::BoundaryID> get_bc_ids() const;
    std::tr1::shared_ptr<libMesh::FunctionBase<Number> > get_func() const;

  private:

    std::vector<GRINS::VariableName> _var_names;
    std::set<GRINS::BoundaryID> _bc_ids;
    std::tr1::shared_ptr<libMesh::FunctionBase<Number> > _func;
    
  };
}
#endif //DBC_CONTAINER_H
