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

#ifndef GRINS_STRESS_STRAIN_LAW_H
#define GRINS_STRESS_STRAIN_LAW_H

// libMesh
#include "libmesh/libmesh_common.h"

// Forward declarations
namespace libMesh
{
  template <typename T>
  class TensorValue;
}

namespace GRINS
{
  class ElasticityTensor;

  template <typename Law>
  class StressStrainLaw
  {
  public:

    StressStrainLaw(){};
    virtual ~StressStrainLaw(){};

    void compute_stress( unsigned int dim,
                         const libMesh::TensorValue<libMesh::Real>& g_contra,
                         const libMesh::TensorValue<libMesh::Real>& g_cov,
                         const libMesh::TensorValue<libMesh::Real>& G_contra,
                         const libMesh::TensorValue<libMesh::Real>& G_cov,
                         libMesh::TensorValue<libMesh::Real>& stress );

    void compute_stress_and_elasticity( unsigned int dim,
                                        const libMesh::TensorValue<libMesh::Real>& g_contra,
                                        const libMesh::TensorValue<libMesh::Real>& g_cov,
                                        const libMesh::TensorValue<libMesh::Real>& G_contra,
                                        const libMesh::TensorValue<libMesh::Real>& G_cov,
                                        libMesh::TensorValue<libMesh::Real>& stress,
                                        ElasticityTensor& C );

    //! This is primarily a helper function for the plane stress cases
    libMesh::Real compute_33_stress( const libMesh::TensorValue<libMesh::Real>& g_contra,
                                     const libMesh::TensorValue<libMesh::Real>& g_cov,
                                     const libMesh::TensorValue<libMesh::Real>& G_contra,
                                     const libMesh::TensorValue<libMesh::Real>& G_cov );
  };

  template <typename Law>
  inline
  void StressStrainLaw<Law>::compute_stress( unsigned int dim,
                                             const libMesh::TensorValue<libMesh::Real>& g_contra,
                                             const libMesh::TensorValue<libMesh::Real>& g_cov,
                                             const libMesh::TensorValue<libMesh::Real>& G_contra,
                                             const libMesh::TensorValue<libMesh::Real>& G_cov,
                                             libMesh::TensorValue<libMesh::Real>& stress )
  {
    static_cast<Law*>(this)->compute_stress_imp(dim,g_contra,g_cov,G_contra,G_cov,stress);
    return;
  }

  template <typename Law>
  inline
  void StressStrainLaw<Law>::compute_stress_and_elasticity( unsigned int dim,
                                                            const libMesh::TensorValue<libMesh::Real>& g_contra,
                                                            const libMesh::TensorValue<libMesh::Real>& g_cov,
                                                            const libMesh::TensorValue<libMesh::Real>& G_contra,
                                                            const libMesh::TensorValue<libMesh::Real>& G_cov,
                                                            libMesh::TensorValue<libMesh::Real>& stress,
                                                            ElasticityTensor& C )
  {
    static_cast<Law*>(this)->compute_stress_and_elasticity_imp(dim,g_contra,g_cov,G_contra,G_cov,stress,C);
    return;
  }

  template <typename Law>
  inline
  libMesh::Real StressStrainLaw<Law>::compute_33_stress( const libMesh::TensorValue<libMesh::Real>& g_contra,
                                                         const libMesh::TensorValue<libMesh::Real>& g_cov,
                                                         const libMesh::TensorValue<libMesh::Real>& G_contra,
                                                         const libMesh::TensorValue<libMesh::Real>& G_cov )
  {
    return static_cast<Law*>(this)->compute_33_stress_imp( g_contra, g_cov, G_contra, G_cov );
  }

} // end namespace GRINS

#endif //GRINS_STRESS_STRAIN_LAW_H
