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

#ifndef GRINS_BOUNDARY_CONDITION_FACTORY_INITIALIZER_H
#define GRINS_BOUNDARY_CONDITION_FACTORY_INITIALIZER_H

namespace GRINS
{
  //! Initialize all Factory objects related to boundary conditions
  /*! To avoid symbol stripping from static linking, we use this
    class to initialize/register the boundary condition factory objects.

    Relevant discussion: http://stackoverflow.com/questions/5202142/static-variable-initialization-over-a-library*/
  class BoundaryConditionFactoryInitializer
  {
  public:
    BoundaryConditionFactoryInitializer();
    ~BoundaryConditionFactoryInitializer(){}
  };
}

#endif // GRINS_BOUNDARY_CONDITION_FACTORY_INITIALIZER_H
