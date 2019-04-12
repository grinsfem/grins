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
#ifndef GRINS_GENERIC_IC_HANDLER_H
#define GRINS_GENERIC_IC_HANDLER_H

//GRINS
#include "grins/variable_name_defaults.h"
#include "grins/var_typedefs.h"
#include "grins/cached_values.h"
#include "grins/ic_handling_base.h"

//libMesh
#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"
#include "libmesh/point.h"
#include "libmesh/function_base.h"

namespace GRINS
{

  //! Base class for reading and handling initial conditions for physics classes
  class GenericICHandler : public ICHandlingBase
  {
  public:

    GenericICHandler(const std::string& physics_name,
                     const GetPot& input);

    virtual ~GenericICHandler();
  };

}
#endif // GRINS_GENERIC_IC_HANDLER
