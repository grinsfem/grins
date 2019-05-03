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

#ifndef GRINS_CARTESIAN_SOLID_MECHANICS_H
#define GRINS_CARTESIAN_SOLID_MECHANICS_H

//GRINS
#include "grins/solid_mechanics_abstract.h"
#include "grins/assembly_context.h"

namespace GRINS
{
  template<unsigned int Dim>
  class CartesianSolidMechanics : public SolidMechanicsAbstract<Dim>
  {
  public:

    CartesianSolidMechanics( const PhysicsName & physics_name,
                             const PhysicsName & core_physics_name,
                             const GetPot & input );

    CartesianSolidMechanics() = delete;

    virtual ~CartesianSolidMechanics() = default;

    //! Initialize context for added physics variables
    virtual void init_context( AssemblyContext & context ) override;

    virtual void mass_residual( bool compute_jacobian, AssemblyContext & context ) override;

  protected:

    //! 2D deformation gradient for plane strain
    /* F(2,2) = 1 will be consistent with plane strain for everything that consumes F */
    libMesh::Tensor form_def_gradient( const libMesh::Gradient & grad_u,
                                       const libMesh::Gradient & grad_v ) const;

    //! 3D deformation gradient
    libMesh::Tensor form_def_gradient( const libMesh::Gradient & grad_u,
                                       const libMesh::Gradient & grad_v,
                                       const libMesh::Gradient & grad_w ) const;

  };

} // end namespace GRINS

#endif // GRINS_CARTESIAN_SOLID_MECHANICS_H
