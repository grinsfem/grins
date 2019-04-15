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

    HyperelasticStrainEnergy() = default;

    virtual ~HyperelasticStrainEnergy() = default;

    //! dW/dI1
    libMesh::Number dI1( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const;

    //! dW/dI2
    libMesh::Number dI2( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const;

    //! dW/dI3
    libMesh::Number dI3( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const;

    //! d^2W/dI1^2
    libMesh::Number dI12( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const;

    //! d^2W/dI2^2
    libMesh::Number dI22( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const;

    //! d^2W/dI3^2
    libMesh::Number dI32( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const;

    //! d^2W/dI1dI2
    libMesh::Number dI1dI2( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const;

    //! d^2W/dI1dI3
    libMesh::Number dI1dI3( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const;

    //! d^2W/dI2dI3
    libMesh::Number dI2dI3( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const;
  };

  template<typename Function>
  libMesh::Number HyperelasticStrainEnergy<Function>::dI1( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const
  {
    return static_cast<const Function*>(this)->dI1_imp(I1,I2,I3);
  }

  template<typename Function>
  libMesh::Number HyperelasticStrainEnergy<Function>::dI2( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const
  {
    return static_cast<const Function*>(this)->dI2_imp(I1,I2,I3);
  }

  template<typename Function>
  libMesh::Number HyperelasticStrainEnergy<Function>::dI3( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const
  {
    return static_cast<const Function*>(this)->dI3_imp(I1,I2,I3);
  }

  template<typename Function>
  libMesh::Number HyperelasticStrainEnergy<Function>::dI12( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const
  {
    return static_cast<const Function*>(this)->dI12_imp(I1,I2,I3);
  }

  template<typename Function>
  libMesh::Number HyperelasticStrainEnergy<Function>::dI22( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const
  {
    return static_cast<const Function*>(this)->dI22_imp(I1,I2,I3);
  }

  template<typename Function>
  libMesh::Number HyperelasticStrainEnergy<Function>::dI32( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const
  {
    return static_cast<const Function*>(this)->dI32_imp(I1,I2,I3);
  }

  template<typename Function>
  libMesh::Number HyperelasticStrainEnergy<Function>::dI1dI2( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const
  {
    return static_cast<const Function*>(this)->dI1dI2_imp(I1,I2,I3);
  }

  template<typename Function>
  libMesh::Number HyperelasticStrainEnergy<Function>::dI1dI3( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const
  {
    return static_cast<const Function*>(this)->dI1dI3_imp(I1,I2,I3);
  }

  template<typename Function>
  libMesh::Number HyperelasticStrainEnergy<Function>::dI2dI3( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const
  {
    return static_cast<const Function*>(this)->dI2dI3_imp(I1,I2,I3);
  }

} // end namespace GRINS

#endif // GRINS_HYPERELASTIC_STRAIN_ENERGY_H
