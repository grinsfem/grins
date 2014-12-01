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
#include "grins/hyperelasticity.h"

// libMesh
#include "libmesh/tensor_value.h"

namespace GRINS
{

  template <typename StrainEnergy>
  Hyperelasticity<StrainEnergy>::Hyperelasticity( const GetPot& input )
    : _W(input)
  {
    return;
  }

  template <typename StrainEnergy>
  Hyperelasticity<StrainEnergy>::~Hyperelasticity()
  {
    return;
  }
  
  template <typename StrainEnergy>
  void Hyperelasticity<StrainEnergy>::compute_stress_imp( unsigned int dim,
                                                          const libMesh::TensorValue<libMesh::Real>& g_contra,
                                                          const libMesh::TensorValue<libMesh::Real>& g_cov,
                                                          const libMesh::TensorValue<libMesh::Real>& G_contra,
                                                          const libMesh::TensorValue<libMesh::Real>& G_cov,
                                                          libMesh::TensorValue<libMesh::Real>& stress )
  {
    stress.zero();

    // Compute strain invariants
    libMesh::Real I3 = (g_contra*G_cov).det();

    libMesh::Real I1 = 0.0;
    libMesh::Real I2 = 0.0;
    for( unsigned int i = 0; i < dim; i++ )
      {
        for( unsigned int j = 0; j < dim; j++ )
          {
            I1 += g_contra(i,j)*G_cov(i,j);
            I2 += G_contra(i,j)*g_cov(i,j);
          }
      }

    I2 *= I3;

    libMesh::Real dWdI1 = _W.dI1(I1,I2,I3);
    libMesh::Real dWdI2 = _W.dI2(I1,I2,I3);
    libMesh::Real dWdI3 = _W.dI3(I1,I2,I3);

    // Now compute stress
    for( unsigned int i = 0; i < dim; i++ )
      {
        for( unsigned int j = 0; j < dim; j++ )
          {
            for( unsigned int k = 0; k < dim; k++ )
              {
                for( unsigned int l = 0; l < dim; l++ )
                  {
                    stress(i,j) += 2.0*dWdI1*g_contra(i,j)
                      + 2.0*dWdI2*(I1*g_contra(i,j) - g_contra(i,k)*g_contra(j,l)*G_cov(k,l))
                      + 2.0*dWdI3*G_contra(i,j);
                  }
              }
          }
      }

    return;
  }

  template <typename StrainEnergy>
  void Hyperelasticity<StrainEnergy>::compute_stress_and_elasticity_imp( unsigned int dim,
                                                                         const libMesh::TensorValue<libMesh::Real>& g_contra,
                                                                         const libMesh::TensorValue<libMesh::Real>& g_cov,
                                                                         const libMesh::TensorValue<libMesh::Real>& G_contra,
                                                                         const libMesh::TensorValue<libMesh::Real>& G_cov,
                                                                         libMesh::TensorValue<libMesh::Real>& stress,
                                                                         ElasticityTensor& C )
  {
    libmesh_not_implemented();
    return;
  }

  template <typename StrainEnergy>
  libMesh::Real Hyperelasticity<StrainEnergy>::compute_33_stress_imp( const libMesh::TensorValue<libMesh::Real>& g_contra,
                                                                      const libMesh::TensorValue<libMesh::Real>& g_cov,
                                                                      const libMesh::TensorValue<libMesh::Real>& G_contra,
                                                                      const libMesh::TensorValue<libMesh::Real>& G_cov )
  {
    libmesh_not_implemented();
    return 0.0;
  }

} // end namespace GRINS
