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

#ifndef GRINS_HOOKES_LAW_1D_H
#define GRINS_HOOKES_LAW_1D_H

// GRINS
#include "grins/parameter_user.h"
#include "grins/stress_strain_law.h"

// Forward declarations
class GetPot;

namespace GRINS
{

  //! Hooke's law specialized for one-dimensional problems
  /*!
   * General form is not conducive to one-dimensional problems
   * so this should be used in those cases.
   */
  class HookesLaw1D : public StressStrainLaw<HookesLaw1D>,
                      public ParameterUser
  {
  public:

    HookesLaw1D( const GetPot& input );

    HookesLaw1D( const GetPot& input, const std::string& material );

    virtual ~HookesLaw1D();

    // So we can make implementations private
    friend class StressStrainLaw<HookesLaw1D>;

  private:

    HookesLaw1D();

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

    //! Lam\'{e} constant
    libMesh::Real _E;

    //! Lam\'{e} constant
    libMesh::Real _nu;

  };

} // end namespace GRINS

#endif // GRINS_HOOKES_LAW_1D_H
