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

#ifndef GRINS_CARTESIAN_STRESS_STRAIN_LAW_H
#define GRINS_CARTESIAN_STRESS_STRAIN_LAW_H

// GRINS
#include "grins/cartesian_nonlinear_mechanics.h"

namespace GRINS
{
  //! Class defining interface to evaluate basic quantities needed for nonlinear elasticity
  /*!
   * At the most general level, we only need the PK1 stress and it's corresponding
   * elasticity tensor, so this minimal class defines the interface for this.
   * It is intended that subclasses take the deformation gradient F at construction time
   * and initialize everything that will be needed to compute the PK1 stress and the
   * componenets of the elasticity tensor. The thinking here is that this object can be
   * constructed within a quadrature loop (when F is computed) and the reused for all
   * computations for that quadrature point.
   */
  template<typename Model>
  class CartesianStressStrainLaw : public CartesianNonlinearMechanics
  {
  public:

    CartesianStressStrainLaw(const libMesh::Tensor & F)
      : _F(F)
    {}

    virtual ~CartesianStressStrainLaw() = default;

    //! Returns first Piola-Kirchoff stress.
    libMesh::Tensor pk1_stress() const;

    //! Returns i,j,k,l entry of elasticity tensor
    libMesh::Number elasticity_tensor( unsigned int i, unsigned int j, unsigned int k, unsigned int l ) const;

  protected:

    //! Cached deformation gradient
    const libMesh::Tensor & _F;
  };

  template <typename Model>
  inline
  libMesh::Tensor CartesianStressStrainLaw<Model>::pk1_stress() const
  {
    return static_cast<const Model*>(this)->pk1_stress_imp();
  }

  template <typename Model>
  inline
  libMesh::Number CartesianStressStrainLaw<Model>::elasticity_tensor
  ( unsigned int i, unsigned int j, unsigned int k, unsigned int l ) const
  {
    return static_cast<const Model*>(this)->elasticity_tensor_imp(i,j,k,l);
  }

} // end namespace GRINS

#endif // GRINS_CARTESIAN_STRESS_STRAIN_LAW_H
