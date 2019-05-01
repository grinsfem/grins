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

#ifndef GRINS_CARTESIAN_NONLINEAR_MECHANICS_H
#define GRINS_CARTESIAN_NONLINEAR_MECHANICS_H

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"

namespace GRINS
{
  //! Helper class for doing computations for nonlinear mechanics
  class CartesianNonlinearMechanics
  {
  public:
    CartesianNonlinearMechanics() = default;
    virtual ~CartesianNonlinearMechanics() = default;

    libMesh::Tensor right_cauchy_green( const libMesh::Tensor & F ) const
    { return F.transpose()*F; }

    void compute_invariants( const libMesh::Tensor & C,
                             libMesh::Number & I1, libMesh::Number & I2, libMesh::Number & I3 ) const;

    libMesh::Real delta( int i, int j ) const
    { return (i==j) ? 1.0 : 0.0; }
  };

  inline
  void CartesianNonlinearMechanics::compute_invariants( const libMesh::Tensor & C,
                                                        libMesh::Number & I1,
                                                        libMesh::Number & I2,
                                                        libMesh::Number & I3 ) const
  {
    I1 = C.tr();

    // I2 = 0.5*( (tr(C))^2 - tr(C^2) )
    I2 = 0.5*( I1*I1 - (C*C).tr() );

    I3 = C.det();
  }

} // end namespace GRINS

#endif // GRINS_CARTESIAN_NONLINEAR_MECHANICS_H
