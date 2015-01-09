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

// GRINS
#include "grins/elasticity_tensor.h"

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
    // We're treating a_* and A_* as 2x2, but we're cheating to pick up lambda^2
    libMesh::Real lambda_sq = A_cov(2,2);

    libMesh::Real A_over_a = 1.0/lambda_sq; // We're incompressible

    libMesh::Real I1, I2;
    this->compute_I1_I2(a_contra,a_cov,A_contra,A_cov,lambda_sq,A_over_a,I1,I2);

    libMesh::Real a_term, A_term;
    this->compute_stress_terms( lambda_sq, A_over_a, I1, I2, a_term, A_term );

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

    // We're treating a_* and A_* as 2x2, but we're cheating to pick up lambda^2
    libMesh::Real lambda_sq = A_cov(2,2);

    libMesh::Real A_over_a = 1.0/lambda_sq; // We're incompressible

    libMesh::Real I1, I2;
    this->compute_I1_I2(a_contra,a_cov,A_contra,A_cov,lambda_sq,A_over_a,I1,I2);

    libMesh::Real a_term, A_term;
    this->compute_stress_terms( lambda_sq, A_over_a, I1, I2, a_term, A_term );

    libMesh::TensorValue<libMesh::Real> daterm_dstrain, dAterm_dstrain;
    this->compute_stress_deriv_terms( lambda_sq, A_over_a, I1, I2, A_contra, daterm_dstrain, dAterm_dstrain );

    ElasticityTensor dAcontra_dstrain;
    this->compute_Acontra_deriv( A_contra, dAcontra_dstrain );

    // Now compute stress
    stress.zero();
    for( unsigned int alpha = 0; alpha < 2; alpha++ )
      {
        for( unsigned int beta = 0; beta < 2; beta++ )
          {
            stress(alpha,beta) = a_contra(alpha,beta)*a_term + A_contra(alpha,beta)*A_term;

            for( unsigned int lambda = 0; lambda < 2; lambda++ )
              {
                for( unsigned int mu = 0; mu < 2; mu++ )
                  {
                    C(alpha,beta,lambda,mu) = a_contra(alpha,beta)*daterm_dstrain(lambda,mu)
                                            + dAcontra_dstrain(alpha,beta,lambda,mu)*A_term
                                            + A_contra(alpha,beta)*dAterm_dstrain(lambda,mu);
                  }
              }
          }
      }

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

  template <typename StrainEnergy>
  void IncompressiblePlaneStressHyperelasticity<StrainEnergy>::compute_stress_terms( libMesh::Real lambda_sq, libMesh::Real A_over_a,
                                                                                     libMesh::Real I1, libMesh::Real I2,
                                                                                     libMesh::Real& a_term, libMesh::Real& A_term ) const
  {
    libMesh::Real dWdI1 = _W.dI1(I1,I2,1.0); // We're incompressible
    libMesh::Real dWdI2 = _W.dI2(I1,I2,1.0);

    // Notation used in Green/Adkins
    // Comes from enforcing plane stress
    libMesh::Real p = -2.0*lambda_sq*( dWdI1 + dWdI2*(I1-lambda_sq) );

    a_term = 2.0*(dWdI1 + dWdI2*lambda_sq);
    A_term = 2.0*dWdI2*A_over_a + p;

    return;
  }

  template <typename StrainEnergy>
  void IncompressiblePlaneStressHyperelasticity<StrainEnergy>::compute_stress_deriv_terms( libMesh::Real lambda_sq, libMesh::Real A_over_a,
                                                                                           libMesh::Real I1, libMesh::Real I2,
                                                                                           const libMesh::TensorValue<libMesh::Real>& A_contra,
                                                                                           libMesh::TensorValue<libMesh::Real>& daterm_dstrain,
                                                                                           libMesh::TensorValue<libMesh::Real>& dAterm_dstrain) const
  {
    daterm_dstrain.zero();
    dAterm_dstrain.zero();

    libMesh::Real dWdI1 = _W.dI1(I1,I2,1.0); // We're incompressible
    libMesh::Real dWdI2 = _W.dI2(I1,I2,1.0);

    // A = det(A_cov) = 1/det(A_contra)
    //libMesh::Real A = 1.0/( A_contra(0,0)*A_contra(1,1) - A_contra(0,1)*A_contra(1,0) );

    for( unsigned int alpha = 0; alpha < 2; alpha++ )
      {
        for( unsigned int beta = 0; beta < 2; beta++ )
          {
            // a_term = 2.0*(dWdI1 + dWdI2*lambda^2);
            // => da_dstrain = 2.0*dWdI2*dlambda^2_dstrain
            // dlambda_sq_dstrain = -2*lambda^2 A^{alpha,beta}
            daterm_dstrain(alpha,beta) = 2.0*dWdI2*(-2.0*lambda_sq*A_contra(alpha,beta));

            // A_term = 2.0*dWdI2*A_over_a + p;
            // A = det(A_cov) ==> dA_dstrain = 2*A*A_contra(alpha,beta)
            // p = -2.0*lambda_sq*( dWdI1 + dWdI2*(I1-lambda_sq) );

            dAterm_dstrain(alpha,beta) = 2.0*dWdI2*A_over_a*(2.0*A_contra(alpha,beta))
                                       - 2.0*(-2.0*lambda_sq*A_contra(alpha,beta))*( dWdI1 + dWdI2*(I1-lambda_sq) )
                                       - 2.0*lambda_sq*( dWdI2*(2.0*lambda_sq*A_contra(alpha,beta)) );
          }
      }

    return;
  }

  template <typename StrainEnergy>
  void IncompressiblePlaneStressHyperelasticity<StrainEnergy>::compute_Acontra_deriv( const libMesh::TensorValue<libMesh::Real>& A_contra,
                                                                                      ElasticityTensor& dAcontra_dstrain ) const
  {
    for( unsigned int alpha = 0; alpha < 2; alpha++ )
      {
        for( unsigned int beta = 0; beta < 2; beta++ )
          {
            for( unsigned int lambda = 0; lambda < 2; lambda++ )
              {
                for( unsigned int mu = 0; mu < 2; mu++ )
                  {
                    dAcontra_dstrain(alpha,beta,lambda,mu) = -0.5*( A_contra(alpha,lambda)*A_contra(beta,mu)
                                                                    + A_contra(alpha,mu)*A_contra(beta,lambda) )
                                                             -A_contra(alpha,beta)*A_contra(lambda,mu);
                  }
              }
          }
      }

    return;
  }

  template <typename StrainEnergy>
  libMesh::Real IncompressiblePlaneStressHyperelasticity<StrainEnergy>::compute_33_stress_imp( const libMesh::TensorValue<libMesh::Real>& /*g_contra*/,
                                                                                               const libMesh::TensorValue<libMesh::Real>& /*g_cov*/,
                                                                                               const libMesh::TensorValue<libMesh::Real>& /*G_contra*/,
                                                                                               const libMesh::TensorValue<libMesh::Real>& /*G_cov*/ )
  {
    std::cerr << "Error: compute_33_stress shouldn't be called for incompressible materials." << std::endl;
    libmesh_error();
    return 0.0;
  }

} // end namespace GRINS
