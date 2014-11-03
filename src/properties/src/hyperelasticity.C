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
  void Hyperelasticity<StrainEnergy>::compute_stress_imp( const libMesh::TensorValue<libMesh::Real>& g_contra,
                                                          const libMesh::TensorValue<libMesh::Real>& G_contra,
                                                          const libMesh::TensorValue<libMesh::Real>& G_cov,
                                                          const libMesh::TensorValue<libMesh::Real>& strain,
                                                          unsigned int dim,
                                                          libMesh::TensorValue<libMesh::Real>& stress )
  {
    stress.zero();

    for( unsigned int i = 0; i < dim; i++ )
      {
        for( unsigned int j = 0; j < dim; j++ )
          {
            for( unsigned int k = 0; k < dim; k++ )
              {
                for( unsigned int l = 0; l < dim; l++ )
                  {
                    libMesh::Real I1 = g_contra(k,l)*G_cov(k,l);
                    libMesh::Real I3 = G.det()/g.det();
                    libMesh::Real I2 = I3*G_contra(k,l)*g_cov(k,l);

                    stress(i,j) += 2.0*_W.dI1(I1,I2,I3)*g_contra(i,j)
                      + 2.0*_W.dI2(I1,I2,I3)*(I1*g_contra(i,j) - g_contra(i,k)*g_contra(j,l)*G_cov(k,l))
                      + 2.0*I3*G_contra(i,j);
                  }
              }
          }
      }

    return;
  }

} // end namespace GRINS
