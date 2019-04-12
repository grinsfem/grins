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

#ifndef GRINS_ONED_CURVILINEAR_SOLID_MECHANICS_H
#define GRINS_ONED_CURVILINEAR_SOLID_MECHANICS_H

//GRINS
#include "grins/solid_mechanics_abstract.h"
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/fe_base.h"

namespace GRINS
{
  class OneDCurvilinearSolidMechanics : public SolidMechanicsAbstract<1>
  {
  public:

    OneDCurvilinearSolidMechanics( const PhysicsName& physics_name, const GetPot& input );

    OneDCurvilinearSolidMechanics() = delete;

    virtual ~OneDCurvilinearSolidMechanics() = default;

    //! Initialize context for added physics variables
    virtual void init_context( AssemblyContext& context );

  protected:

    //! Cross-sectional area of the cable
    libMesh::Real _A;

  };

}

#endif // GRINS_ONED_CURVILINEAR_SOLID_MECHANICS_H
