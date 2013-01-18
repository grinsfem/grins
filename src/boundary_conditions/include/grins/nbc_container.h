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

#ifndef GRINS_NBC_CONTAINER_H
#define GRINS_NBC_CONTAINER_H

#include "grins/neumann_func_obj.h"

namespace GRINS
{
  //! Simple helper class to setup general Neumann boundary conditions
  /*! This class is to temporarily stash data necessary for setting
      up GRINS::NeumannFuncObj objects. Actual instantiation
      of GRINS::NeumannFuncObj object is handled internally by
      BCHandling objects. This class is structured on a per-bc ID basis.
      That is, each object works handles functions for one boundary ID at a time,
      but can handle multiple variable/function pairs. Currently, it is assumed
      one function per variable. */
  class NBCContainer
  {
  public:
    NBCContainer();
    ~NBCContainer();

    //! Add variable for which this boundary condition is to be applied. 
    void set_bc_id( BoundaryID bc_id );

    //! Add boundary id and corresponding functor object to be applied on that boundary
    void add_var_func_pair( VariableIndex var, 
			    std::tr1::shared_ptr<NeumannFuncObj> func );

    BoundaryID get_bc_id() const;

    std::tr1::shared_ptr<NeumannFuncObj> get_func( VariableIndex var ) const;

  protected:
    
    BoundaryID _bc_id;
    std::map<VariableIndex,std::tr1::shared_ptr<NeumannFuncObj> > _funcs;

  };

} // namespace GRINS

#endif // GRINS_NBC_CONTAINER_H
