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

#ifndef GRINS_HOOKES_LAW_H
#define GRINS_HOOKES_LAW_H

// GRINS
#include "grins/stress_strain_law.h"
#include "grins/elasticity_tensor.h"

// Forward declarations
class GetPot;

namespace GRINS
{

  //! Elasticity tensor for Hooke's law
  /*!
   * Uses Lam\'{e} constants, but can parse Young's modulus and Poisson's ratio
   * if desired. By default, is constructed for Cartesian coordinate systems. If
   * working with curvilinear coordinate systems, the user should call the
   * set_deformation method before calling operator().
   */
  class HookesLaw : public StressStrainLaw<HookesLaw>
  {
  public:

    HookesLaw( const GetPot& input );
    virtual ~HookesLaw();

    // So we can make implementations private
    friend class StressStrainLaw<HookesLaw>;

  private:

    HookesLaw();

    //! Parse properties from input
    void read_input_options(const GetPot& input);

    void compute_stress_imp( unsigned int dim,
                             const libMesh::TensorValue<libMesh::Real>& g_contra,
                             const libMesh::TensorValue<libMesh::Real>& g_cov,
                             const libMesh::TensorValue<libMesh::Real>& G_contra,
                             const libMesh::TensorValue<libMesh::Real>& G_cov,
                             libMesh::TensorValue<libMesh::Real>& stress );

    void compute_stress_and_elasticity_imp( unsigned int dim,
                                            const libMesh::TensorValue<libMesh::Real>& g_contra,
                                            const libMesh::TensorValue<libMesh::Real>& g_cov,
                                            const libMesh::TensorValue<libMesh::Real>& G_contra,
                                            const libMesh::TensorValue<libMesh::Real>& G_cov,
                                            libMesh::TensorValue<libMesh::Real>& stress,
                                            ElasticityTensor& C );

    libMesh::Real compute_33_stress_imp( const libMesh::TensorValue<libMesh::Real>& g_contra,
                                         const libMesh::TensorValue<libMesh::Real>& g_cov,
                                         const libMesh::TensorValue<libMesh::Real>& G_contra,
                                         const libMesh::TensorValue<libMesh::Real>& G_cov );

    ElasticityTensor _C;

    //! Lam\'{e} constant
    libMesh::Real _lambda;

    //! Lam\'{e} constant
    libMesh::Real _mu;

  };

} // end namespace GRINS

#endif // GRINS_HOOKES_LAW_H
