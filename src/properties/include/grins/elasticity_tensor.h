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

#ifndef GRINS_ELASTICITY_TENSOR_H
#define GRINS_ELASTICITY_TENSOR_H

// libMesh
#include "libmesh/libmesh_common.h"

namespace GRINS
{
  class ElasticityTensor
  {
  public:

    ElasticityTensor(){};
    virtual ~ElasticityTensor(){};

    //! Value of C_{ijkl}
    libMesh::Real operator()( unsigned int i, unsigned int j, unsigned int k, unsigned int l ) const;

    libMesh::Real& operator()( unsigned int i, unsigned int j, unsigned int k, unsigned int l );

  protected:

    //! Elasticity tensor
    /*! \todo Does this guarantee continuous memory? If not, we should use a datastructure that does */
    libMesh::Real _C[3][3][3][3];
  };

  inline
  libMesh::Real ElasticityTensor::operator()( unsigned int i, unsigned int j, unsigned int k, unsigned int l ) const
  {
    return _C[i][j][k][l];
  }

  inline
  libMesh::Real& ElasticityTensor::operator()( unsigned int i, unsigned int j, unsigned int k, unsigned int l )
  {
    return _C[i][j][k][l];
  }

} // end namespace GRINS

#endif // GRINS_ELASTICITY_TENSOR_H
