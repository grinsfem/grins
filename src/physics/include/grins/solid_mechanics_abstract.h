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

#ifndef GRINS_SOLID_MECHANICS_ABSTRACT_H
#define GRINS_SOLID_MECHANICS_ABSTRACT_H

//GRINS
#include "grins/physics.h"
#include "grins/multi_component_vector_variable.h"

// libMesh
#include "libmesh/fem_context.h"

namespace GRINS
{
  template<unsigned int Dim>
  class SolidMechanicsAbstract : public Physics
  {
  public:

    SolidMechanicsAbstract( const PhysicsName & physics_name,
                            const PhysicsName & core_physics_name,
                            const GetPot & input );

    virtual ~SolidMechanicsAbstract() =default;

    virtual void set_time_evolving_vars( libMesh::FEMSystem* system ) override;

  protected:

    DisplacementVariable& _disp_vars;

    //! Solid density
    libMesh::Real _rho;

    const libMesh::FEGenericBase<libMesh::Real> * get_fe( const AssemblyContext & context );

    typedef const libMesh::DenseSubVector<libMesh::Number>& (libMesh::DiffContext::*VarFuncType)(unsigned int) const;

    typedef void (libMesh::FEMContext::*InteriorFuncType)(unsigned int, unsigned int, libMesh::Real&) const;

    typedef libMesh::Real (libMesh::DiffContext::*VarDerivType)() const;

  };

  template<unsigned int Dim>
  inline
  const libMesh::FEGenericBase<libMesh::Real>* SolidMechanicsAbstract<Dim>::get_fe( const AssemblyContext & context )
  {
    // For this Physics, we need to make sure that we grab only the 1D elements
    return context.get_element_fe(_disp_vars.u(),Dim);
  }

} // end namespace GRINS

#endif // GRINS_SOLID_MECHANICS_ABSTRACT_H
