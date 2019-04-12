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

#ifndef GRINS_HYPERELASTIC_STRAIN_ENERGY_H
#define GRINS_HYPERELASTIC_STRAIN_ENERGY_H

// libMesh
#include "libmesh/libmesh_common.h"

namespace GRINS
{
  template<typename Function>
  class HyperelasticStrainEnergy
  {
  public:
    HyperelasticStrainEnergy(){};
    virtual ~HyperelasticStrainEnergy(){};

    libMesh::Real dI1( libMesh::Real I1, libMesh::Real I2, libMesh::Real I3 ) const;
    libMesh::Real dI2( libMesh::Real I1, libMesh::Real I2, libMesh::Real I3 ) const;
    libMesh::Real dI3( libMesh::Real I1, libMesh::Real I2, libMesh::Real I3 ) const;

  };

  template<typename Function>
  libMesh::Real HyperelasticStrainEnergy<Function>::dI1( libMesh::Real I1, libMesh::Real I2, libMesh::Real I3 ) const
  {
    return static_cast<const Function*>(this)->dI1_imp(I1,I2,I3);
  }

  template<typename Function>
  libMesh::Real HyperelasticStrainEnergy<Function>::dI2( libMesh::Real I1, libMesh::Real I2, libMesh::Real I3 ) const
  {
    return static_cast<const Function*>(this)->dI2_imp(I1,I2,I3);
  }

  template<typename Function>
  libMesh::Real HyperelasticStrainEnergy<Function>::dI3( libMesh::Real I1, libMesh::Real I2, libMesh::Real I3 ) const
  {
    return static_cast<const Function*>(this)->dI3_imp(I1,I2,I3);
  }

} // end namespace GRINS



#endif // GRINS_HYPERELASTIC_STRAIN_ENERGY_H
