//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_SOLID_MECHANICS_ABSTRACT_H
#define GRINS_SOLID_MECHANICS_ABSTRACT_H

//GRINS
#include "grins/physics.h"
#include "grins/multi_component_vector_variable.h"

// libMesh
#include "libmesh/fem_context.h"

namespace GRINS
{
  class SolidMechanicsAbstract : public Physics
  {
  public:

    SolidMechanicsAbstract( const GRINS::PhysicsName& physics_name,
                            const GetPot& input );

    virtual ~SolidMechanicsAbstract(){};

    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

  protected:

    DisplacementVariable& _disp_vars;

    typedef const libMesh::DenseSubVector<libMesh::Number>& (libMesh::DiffContext::*VarFuncType)(unsigned int) const;

    typedef void (libMesh::FEMContext::*InteriorFuncType)(unsigned int, unsigned int, libMesh::Real&) const;

    typedef libMesh::Real (libMesh::DiffContext::*VarDerivType)() const;

  private:

    SolidMechanicsAbstract();

  };

} // end namespace GRINS

#endif // GRINS_SOLID_MECHANICS_ABSTRACT_H
