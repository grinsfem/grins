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

#ifndef GRINS_THREED_SOLID_MECHANICS_BASE_H
#define GRINS_THREED_SOLID_MECHANICS_BASE_H

//GRINS
#include "grins/solid_mechanics_abstract.h"
#include "grins/assembly_context.h"

namespace GRINS
{
  class ThreeDSolidMechanicsBase : public SolidMechanicsAbstract<3>
  {
  public:

    ThreeDSolidMechanicsBase( const PhysicsName & physics_name,
                              const PhysicsName & core_physics_name,
                              const GetPot & input );

    ThreeDSolidMechanicsBase() = delete;

    virtual ~ThreeDSolidMechanicsBase() = default;

    //! Initialize context for added physics variables
    virtual void init_context( AssemblyContext & context );

    virtual void mass_residual( bool compute_jacobian, AssemblyContext & context ) override;

  protected:

    libMesh::Tensor form_def_gradient( const libMesh::Gradient & grad_u,
                                       const libMesh::Gradient & grad_v,
                                       const libMesh::Gradient & grad_w ) const;

    libMesh::Tensor compute_right_cauchy_def( const libMesh::Tensor & F ) const
    { return F.transpose()*F; }

    void compute_invariants( const libMesh::Tensor & C,
                             libMesh::Number & I1, libMesh::Number & I2, libMesh::Number & I3 ) const;

    libMesh::Real delta( int i, int j ) const
    { return (i==j) ? 1.0 : 0.0; }

  };

} // end namespace GRINS

#endif // GRINS_THREED_SOLID_MECHANICS_BASE_H
