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

#ifndef GRINS_COMPRESSIBLE_MOONEY_RIVLIN_H
#define GRINS_COMPRESSIBLE_MOONEY_RIVLIN_H

// GRINS
#include "grins/hyperelastic_strain_energy.h"
#include "grins/parameter_user.h"

// Forward declarations
class GetPot;

namespace GRINS
{
  class CompressibleMooneyRivlin : public HyperelasticStrainEnergy<CompressibleMooneyRivlin>,
                                   public ParameterUser
  {
  public:
    CompressibleMooneyRivlin( const GetPot & input );

    CompressibleMooneyRivlin() = delete;

    CompressibleMooneyRivlin( const GetPot & input, const std::string & material );

    virtual ~CompressibleMooneyRivlin() = default;

    // So we can make implementation private
    friend class  HyperelasticStrainEnergy<CompressibleMooneyRivlin>;

  private:

    libMesh::Number dI1_imp( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const;
    libMesh::Number dI2_imp( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const;
    libMesh::Number dI3_imp( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const;

    libMesh::Number dI12_imp( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const;
    libMesh::Number dI22_imp( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const;
    libMesh::Number dI32_imp( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const;

    libMesh::Real _C1;
    libMesh::Real _C2;
    libMesh::Real _C3;
    libMesh::Number dI1dI2_imp( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const;
    libMesh::Number dI1dI3_imp( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const;
    libMesh::Number dI2dI3_imp( libMesh::Number I1, libMesh::Number I2, libMesh::Number I3 ) const;

  };

  inline
  libMesh::Number CompressibleMooneyRivlin::dI1_imp
  ( libMesh::Number /*I1*/, libMesh::Number /*I2*/, libMesh::Number /*I3*/ ) const
  {
    return _C1;
  }

  inline
  libMesh::Number CompressibleMooneyRivlin::dI2_imp
  ( libMesh::Number /*I1*/, libMesh::Number /*I2*/, libMesh::Number /*I3*/ ) const
  {
    return _C2;
  }

  inline
  libMesh::Number CompressibleMooneyRivlin::dI3_imp
  ( libMesh::Number /*I1*/, libMesh::Number /*I2*/, libMesh::Number I3 ) const
  {
    return -(_C1 + 2*_C2 + _C3*std::sqrt(I3) - _C3*I3)/I3;
  }

  inline
  libMesh::Number CompressibleMooneyRivlin::dI12_imp
  ( libMesh::Number /*I1*/, libMesh::Number /*I2*/, libMesh::Number /*I3*/ ) const
  {
    return 0.0;
  }

  inline
  libMesh::Number CompressibleMooneyRivlin::dI22_imp
  ( libMesh::Number /*I1*/, libMesh::Number /*I2*/, libMesh::Number /*I3*/ ) const
  {
    return 0.0;
  }

  inline
  libMesh::Number CompressibleMooneyRivlin::dI32_imp
  ( libMesh::Number /*I1*/, libMesh::Number /*I2*/, libMesh::Number I3 ) const
  {
    return (2*_C1 + 4*_C2 + _C3*std::sqrt(I3))/(2*I3*I3);
  }

  inline
  libMesh::Number CompressibleMooneyRivlin::dI1dI2_imp
  ( libMesh::Number /*I1*/, libMesh::Number /*I2*/, libMesh::Number /*I3*/ ) const
  {
    return 0.0;
  }

  inline
  libMesh::Number CompressibleMooneyRivlin::dI1dI3_imp
  ( libMesh::Number /*I1*/, libMesh::Number /*I2*/, libMesh::Number /*I3*/ ) const
  {
    return 0.0;
  }

  inline
  libMesh::Number CompressibleMooneyRivlin::dI2dI3_imp
  ( libMesh::Number /*I1*/, libMesh::Number /*I2*/, libMesh::Number /*I3*/ ) const
  {
    return 0.0;
  }

} // end namespace GRINS

#endif // GRINS_COMPRESSIBLE_MOONEY_RIVLIN_H
