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

// This class
#include "grins/compressible_hyperelasticity.h"

// GRINS
#include "grins/materials_parsing.h"
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{
  template<typename StrainEnergy>
  CompressibleHyperelasticity<StrainEnergy>::CompressibleHyperelasticity
  ( const PhysicsName & physics_name, const GetPot & input )
    : ThreeDSolidMechanicsBase(physics_name,PhysicsNaming::compressible_hyperelasticity(),input),
      _strain_energy(nullptr)
  {
    const std::string material =
      MaterialsParsing::material_name(input,PhysicsNaming::compressible_hyperelasticity());

    _strain_energy.reset(new StrainEnergy(input,material));
  }

  template<typename StrainEnergy>
  libMesh::Tensor CompressibleHyperelasticity<StrainEnergy>::compute_pk2_stress
  ( const libMesh::Tensor & C,
    const libMesh::Tensor & Cinv,
    const libMesh::Number & I1,
    const libMesh::Number & I3,
    const libMesh::Number & dWdI1,
    const libMesh::Number & dWdI2,
    const libMesh::Number & dWdI3 ) const
  {
    const int dim = 3;
    libMesh::Tensor S;

    for( int i = 0; i < dim; i++ )
      {
        for (int j = 0; j < dim; j++ )
          {
            libMesh::Real dij = this->delta(i,j);
            S(i,j) = 2*( dWdI1*dij + dWdI2*(I1*dij-C(i,j)) + dWdI3*(I3*Cinv(i,j)) );
          }
      }

    return S;
  }
} // end namespace GRINS
