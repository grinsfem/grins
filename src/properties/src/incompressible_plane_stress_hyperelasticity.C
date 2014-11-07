//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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
#include "grins/incompressible_plane_stress_hyperelasticity.h"

// libMesh
#include "libmesh/tensor_value.h"

namespace GRINS
{
  template <typename StrainEnergy>
  IncompressiblePlaneStressHyperelasticity<StrainEnergy>::IncompressiblePlaneStressHyperelasticity( const GetPot& input )
    : _W(input)
  {
    return;
  }

  template <typename StrainEnergy>
  IncompressiblePlaneStressHyperelasticity<StrainEnergy>::~IncompressiblePlaneStressHyperelasticity()
  {
    return;
  }

  template <typename StrainEnergy>
  void IncompressiblePlaneStressHyperelasticity<StrainEnergy>::compute_stress_imp( unsigned int /*dim*/,
                                                                                   const libMesh::TensorValue<libMesh::Real>& a_contra,
                                                                                   const libMesh::TensorValue<libMesh::Real>& a_cov,
                                                                                   const libMesh::TensorValue<libMesh::Real>& A_contra,
                                                                                   const libMesh::TensorValue<libMesh::Real>& A_cov,
                                                                                   libMesh::TensorValue<libMesh::Real>& stress )
  {
    // We're treating a/A as 2x2, but we're cheating to pick up lambda^2
    libMesh::Real lambda_sq = A_cov(2,2);

    libMesh::Real A_over_a = 1.0/lambda_sq; // We're incompressible

    libMesh::Real I1, I2;
    this->compute_I1_I2(a_contra,a_cov,A_contra,A_cov,lambda_sq,A_over_a,I1,I2);

    libMesh::Real dWdI1 = _W.dI1(I1,I2,1.0); // We're incompressible
    libMesh::Real dWdI2 = _W.dI2(I1,I2,1.0);

    // Notation used in Green/Adkins
    // Comes from enforcing plane stress
    libMesh::Real p = -2.0*lambda_sq*( dWdI1 + dWdI2*(I1-lambda_sq) );

    libMesh::Real a_term = 2.0*(dWdI1 + dWdI2*lambda_sq);
    libMesh::Real A_term = 2.0*dWdI2*A_over_a + p;

    // Now compute stress
    stress.zero();
    for( unsigned int alpha = 0; alpha < 2; alpha++ )
      {
        for( unsigned int beta = 0; beta < 2; beta++ )
          {
            stress(alpha,beta) = a_contra(alpha,beta)*a_term + A_contra(alpha,beta)*A_term;
          }
      }

    return;
  }

  template <typename StrainEnergy>
  void IncompressiblePlaneStressHyperelasticity<StrainEnergy>::compute_stress_and_elasticity_imp( unsigned int /*dim*/,
                                                                                                  const libMesh::TensorValue<libMesh::Real>& a_contra,
                                                                                                  const libMesh::TensorValue<libMesh::Real>& a_cov,
                                                                                                  const libMesh::TensorValue<libMesh::Real>& A_contra,
                                                                                                  const libMesh::TensorValue<libMesh::Real>& A_cov,
                                                                                                  libMesh::TensorValue<libMesh::Real>& stress,
                                                                                                  ElasticityTensor& C)
  {
    libmesh_not_implemented();
    return;
  }

  template <typename StrainEnergy>
  void IncompressiblePlaneStressHyperelasticity<StrainEnergy>::compute_I1_I2( const libMesh::TensorValue<libMesh::Real>& a_contra,
                                                                              const libMesh::TensorValue<libMesh::Real>& a_cov,
                                                                              const libMesh::TensorValue<libMesh::Real>& A_contra,
                                                                              const libMesh::TensorValue<libMesh::Real>& A_cov,
                                                                              libMesh::Real lambda_sq, libMesh::Real A_over_a,
                                                                              libMesh::Real& I1, libMesh::Real& I2 ) const
  {
    I1 = 0.0;
    I2 = 0.0;
    for( unsigned int alpha = 0; alpha < 2; alpha++ )
      {
        for( unsigned int beta = 0; beta < 2; beta++ )
          {
            I1 += a_contra(alpha,beta)*A_cov(alpha,beta);
            I2 += a_cov(alpha,beta)*A_contra(alpha,beta);
          }
      }

    I1 += lambda_sq;
    I2 += A_over_a;

    return;
  }

} // end namespace GRINS
