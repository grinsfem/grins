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

#ifndef GRINS_PLANE_STRESS_HYPERELASTICITY_H
#define GRINS_PLANE_STRESS_HYPERELASTICITY_H

// GRINS
#include "grins/stress_strain_law.h"

// libMesh
#include "libmesh/libmesh_common.h"

// Forward declarations
class GetPot;
namespace libMesh
{
  template <typename T>
  class TensorValue;
}

namespace GRINS
{
  template <typename StrainEnergy>
  class IncompressiblePlaneStressHyperelasticity : public StressStrainLaw<IncompressiblePlaneStressHyperelasticity<StrainEnergy> >
  {
  public:
    IncompressiblePlaneStressHyperelasticity( const GetPot& input );
    virtual ~IncompressiblePlaneStressHyperelasticity();

    // So we can make implementation private
    friend class StressStrainLaw<IncompressiblePlaneStressHyperelasticity<StrainEnergy> >;

  private:

    void compute_stress_imp( unsigned int dim,
                             const libMesh::TensorValue<libMesh::Real>& a_contra,
                             const libMesh::TensorValue<libMesh::Real>& a_cov,
                             const libMesh::TensorValue<libMesh::Real>& A_contra,
                             const libMesh::TensorValue<libMesh::Real>& A_cov,
                             libMesh::TensorValue<libMesh::Real>& stress );

    void compute_stress_and_elasticity_imp( unsigned int dim,
                                            const libMesh::TensorValue<libMesh::Real>& a_contra,
                                            const libMesh::TensorValue<libMesh::Real>& a_cov,
                                            const libMesh::TensorValue<libMesh::Real>& A_contra,
                                            const libMesh::TensorValue<libMesh::Real>& A_cov,
                                            libMesh::TensorValue<libMesh::Real>& stress,
                                            ElasticityTensor& C );

    void compute_I1_I2( const libMesh::TensorValue<libMesh::Real>& a_contra,
                        const libMesh::TensorValue<libMesh::Real>& a_cov,
                        const libMesh::TensorValue<libMesh::Real>& A_contra,
                        const libMesh::TensorValue<libMesh::Real>& A_cov,
                        libMesh::Real lambda_sq, libMesh::Real A_over_a,
                        libMesh::Real& I1, libMesh::Real& I2) const;

    void compute_stress_terms( libMesh::Real lambda_sq, libMesh::Real A_over_a,
                               libMesh::Real I1, libMesh::Real I2,
                               libMesh::Real& a_term, libMesh::Real& A_term ) const;

    void compute_stress_deriv_terms( libMesh::Real lambda_sq, libMesh::Real A_over_a,
                                     libMesh::Real I1, libMesh::Real I2,
                                     const libMesh::TensorValue<libMesh::Real>& A_contra,
                                     libMesh::TensorValue<libMesh::Real>& daterm_dstrain,
                                     libMesh::TensorValue<libMesh::Real>& dAterm_dstrain) const;

    void compute_Acontra_deriv( const libMesh::TensorValue<libMesh::Real>& A_contra,
                                ElasticityTensor& dAcontra_dstrain ) const;

    libMesh::Real compute_33_stress_imp( const libMesh::TensorValue<libMesh::Real>& g_contra,
                                         const libMesh::TensorValue<libMesh::Real>& g_cov,
                                         const libMesh::TensorValue<libMesh::Real>& G_contra,
                                         const libMesh::TensorValue<libMesh::Real>& G_cov );

    StrainEnergy _W;

  };

} // end namespace GRINS

#endif // GRINS_PLANE_STRESS_HYPERELASTICITY_H
