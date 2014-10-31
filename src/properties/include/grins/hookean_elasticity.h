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

#ifndef GRINS_HOOKEAN_ELASTICITY_H
#define GRINS_HOOKEAN_ELASTICITY_H

// GRINS
#include "grins/elasticity_tensor.h"

// Forward declarations
class GetPot;
namespace libMesh
{
  template<typename T>
  class TensorValue;
}

namespace GRINS
{

  //! Elasticity tensor for Hooke's law
  /*!
   * Uses Lam\'{e} constants, but can parse Young's modulus and Poisson's ratio
   * if desired. By default, is constructed for Cartesian coordinate systems. If
   * working with curvilinear coordinate systems, the user should call the
   * set_deformation method before calling operator().
   */
  class HookeanElasiticty : public ElasticityTensor
  {
  public:

    HookeanElasiticty( const GetPot& input );
    virtual ~HookeanElasiticty();

    //! Compute the elasiticity tensor given the deformation
    /*!
     * This is needed for curvilinear coordinates, shell formulations, etc.
     * By default, we initialize using the Kronecker delta so this only
     * needs to be called when Kronecker delta does not apply.
     */
    void recompute_elasticity( libMesh::TensorValue<libMesh::Real>& g );

  private:

    HookeanElasiticty();

    //! Parse properties from input
    void read_input_options(const GetPot& input);

    //! Lam\'{e} constant
    libMesh::Real _lambda;

    //! Lam\'{e} constant
    libMesh::Real _mu;

  };

} // end namespace GRINS

#endif // GRINS_HOOKEAN_ELASTICITY_H
