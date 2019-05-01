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

#ifndef GRINS_COMPRESSIBLE_HYPERELASTICITY_H
#define GRINS_COMPRESSIBLE_HYPERELASTICITY_H

// GRINS
#include "grins/cartesian_solid_mechanics.h"
#include "grins/hyperelastic_strain_energy.h"

namespace GRINS
{
  template<typename StrainEnergy>
  class CompressibleHyperelasticity : public CartesianSolidMechanics
  {
  public:

    CompressibleHyperelasticity( const PhysicsName & physics_name, const GetPot & input );
    CompressibleHyperelasticity() = delete;

    virtual ~CompressibleHyperelasticity() = default;

    virtual void element_time_derivative( bool compute_jacobian, AssemblyContext & context ) override;

  private:

    libMesh::Tensor compute_pk2_stress( const libMesh::Tensor & C, const libMesh::Tensor & Cinv,
                                        const libMesh::Number & I1, const libMesh::Number & I3,
                                        const libMesh::Number & dWdI1, const libMesh::Number & dWdI2,
                                        const libMesh::Number & dWdI3 ) const;

    libMesh::Number elasticity_tensor( int i, int j, int k, int l,
                                       const libMesh::Tensor & C,
                                       const libMesh::Tensor & Cinv,
                                       const libMesh::Number I1,
                                       const libMesh::Number I2,
                                       const libMesh::Number I3,
                                       const libMesh::Number dWdI2,
                                       const libMesh::Number dWdI3 ) const;

    std::unique_ptr<HyperelasticStrainEnergy<StrainEnergy>> _strain_energy;

  };

} // end namespace GRINS

#endif // GRINS_COMPRESSIBLE_HYPERELASTICITY_H
