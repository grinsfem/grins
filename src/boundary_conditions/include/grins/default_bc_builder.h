//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_DEFAULT_BC_BUILDER_H
#define GRINS_DEFAULT_BC_BUILDER_H

// GRINS
#include "grins/bc_builder.h"

namespace GRINS
{
  //! Manages runtime construction of Dirichlet boundary conditions
  /*! This will parse the input for the request Dirichlet boundary
      conditions and manage their construction. Actual construction of
      the DirichletBoundary objects is delegated to factory
      classes. This builder classes merely manages tasks around the
      factories as needed.  To add new Dirichlet boundary conditions,
      the user should instantiate an appropriate factory sub class. */
  class DefaultBCBuilder : public BCBuilder
  {
  public:

    DefaultBCBuilder()
      : BCBuilder()
    {};

    ~DefaultBCBuilder(){};

  protected:

    virtual void build_bcs( const GetPot& input, MultiphysicsSystem& system,
                            std::vector<SharedPtr<NeumannBCContainer> >& neumann_bcs );

  };
} // end namespace GRINS

#endif // GRINS_DEFAULT_BC_BUILDER_H
