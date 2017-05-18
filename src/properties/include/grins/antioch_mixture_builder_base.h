//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_ANTIOCH_MIXTURE_BUILDER_BASE_H
#define GRINS_ANTIOCH_MIXTURE_BUILDER_BASE_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/vector_utils.h"
#include "antioch/chemical_mixture.h"
#include "antioch/reaction_set.h"
#include "antioch/nasa_mixture.h"
#include "antioch/cea_curve_fit.h"

// C++
#include <string>

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  //! Base class building Antioch mixture wrappers
  /*! This class only worries about building the kinetics
      and the thermo associated with kinetics. Subclasses
      will handle thermo and transport. */
  class AntiochMixtureBuilderBase
  {
  public:
    AntiochMixtureBuilderBase(){}
    ~AntiochMixtureBuilderBase(){}

  };
} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_MIXTURE_BUILDER_BASE_H
