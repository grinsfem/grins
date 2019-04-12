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

#ifndef GRINS_MOONEY_RIVLIN_H
#define GRINS_MOONEY_RIVLIN_H

// GRINS
#include "grins/hyperelastic_strain_energy.h"
#include "grins/parameter_user.h"

// Forward declarations
class GetPot;

namespace GRINS
{
  class MooneyRivlin : public HyperelasticStrainEnergy<MooneyRivlin>,
                       public ParameterUser
  {
  public:
    MooneyRivlin( const GetPot& input );

    MooneyRivlin( const GetPot& input, const std::string& material );

    virtual ~MooneyRivlin();

    // So we can make implementation private
    friend class  HyperelasticStrainEnergy<MooneyRivlin>;

  private:

    MooneyRivlin();

    libMesh::Real dI1_imp( libMesh::Real I1, libMesh::Real I2, libMesh::Real I3 ) const;
    libMesh::Real dI2_imp( libMesh::Real I1, libMesh::Real I2, libMesh::Real I3 ) const;
    libMesh::Real dI3_imp( libMesh::Real I1, libMesh::Real I2, libMesh::Real I3 ) const;

    libMesh::Real _C1;
    libMesh::Real _C2;

  };

  inline
  libMesh::Real MooneyRivlin::dI1_imp( libMesh::Real /*I1*/, libMesh::Real /*I2*/, libMesh::Real /*I3*/ ) const
  {
    return _C1;
  }

  inline
  libMesh::Real MooneyRivlin::dI2_imp( libMesh::Real /*I1*/, libMesh::Real /*I2*/, libMesh::Real /*I3*/ ) const
  {
    return _C2;
  }

  inline
  libMesh::Real MooneyRivlin::dI3_imp( libMesh::Real /*I1*/, libMesh::Real /*I2*/, libMesh::Real /*I3*/ ) const
  {
    return 0.0;
  }


} // end namespace GRINS

#endif // GRINS_MOONEY_RIVLIN_H
