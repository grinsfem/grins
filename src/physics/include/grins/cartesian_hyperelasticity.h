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

#ifndef GRINS_CARTESIAN_HYPERELASTICITY_H
#define GRINS_CARTESIAN_HYPERELASTICITY_H

// GRINS
#include "grins/cartesian_stress_strain_law.h"
#include "grins/hyperelastic_strain_energy.h"

namespace GRINS
{
  //! Stress-strain law for Cartesian hyperelasticity
  /*!
   * This class was meant to be used within a quadrature loop. At each quadrature point,
   * it will compute all the necessary ingredients to compute the pk1 stress and the elasticity
   * tensor at construction time. Thus, the user can obtain a reference to the PK2 stress that will
   * persist the life of this class, etc. It's templated on the strain energy function so the user
   * need only instantiate with a new subclass of HyperelasticStrainEnergy<StrainEnergy>.
   */
  template<typename StrainEnergy>
  class CartesianHyperlasticity : public CartesianStressStrainLaw<CartesianHyperlasticity<StrainEnergy>>
  {
  public:

    CartesianHyperlasticity( const libMesh::Tensor & F, const HyperelasticStrainEnergy<StrainEnergy> & W );

    CartesianHyperlasticity() = delete;
    virtual ~CartesianHyperlasticity() = default;

    // So we can make implementation private
    friend class CartesianStressStrainLaw<CartesianHyperlasticity<StrainEnergy>>;

    const libMesh::Tensor & get_pk2_stress() const
    { return _S;}

    const libMesh::Tensor & get_C_inverse() const
    { return _Cinv; }

    libMesh::Number get_J() const
    { return std::sqrt(_I3); }

  private:

    const HyperelasticStrainEnergy<StrainEnergy> & _W;

    libMesh::Tensor _C, _Cinv;
    libMesh::Number _I1, _I2, _I3;
    libMesh::Gradient _dW;
    libMesh::Tensor _d2W;

    libMesh::Tensor _S;

    libMesh::Tensor pk1_stress_imp() const;

    libMesh::Number elasticity_tensor_imp( unsigned int i, unsigned int j, unsigned int k, unsigned int l ) const;

    libMesh::Tensor compute_pk2_stress() const;
  };

  template<typename StrainEnergy>
  inline
  libMesh::Tensor CartesianHyperlasticity<StrainEnergy>::pk1_stress_imp() const
  {
    return (this->_F)*(this->_S);
  }

  template<typename StrainEnergy>
  inline
  libMesh::Number CartesianHyperlasticity<StrainEnergy>::elasticity_tensor_imp
  ( unsigned int i, unsigned int j, unsigned int k, unsigned int l ) const
  {
    const libMesh::Real dij = this->delta(i,j);
    const libMesh::Real dkl = this->delta(k,l);
    const libMesh::Real dik = this->delta(i,k);
    const libMesh::Real djl = this->delta(j,l);
    const libMesh::Real dil = this->delta(i,l);
    const libMesh::Real djk = this->delta(j,k);

    libMesh::Number Cijkl =
      _d2W(0,0)*(dij*dkl) +
      _d2W(1,1)*( (_I1*dij - _C(i,j))*(_I1*dkl - _C(k,l)) ) +
      _d2W(2,2)*_I3*_I3*_Cinv(i,j)*_Cinv(k,l) +
      _d2W(0,1)*( dij*(_I1*dkl - _C(k,l)) + dkl*(_I1*dij - _C(i,j)) ) +
      _d2W(0,2)*_I3*(dij*_Cinv(k,l) + _Cinv(i,j)*dkl) +
      _d2W(1,2)*_I3*( _Cinv(k,l)*(_I1*dij - _C(i,j)) + _Cinv(i,j)*(_I1*dkl-_C(k,l)) ) +
      _dW(1)*(dij*dkl - 0.5*(dik*djl + dil*djk) ) +
      _dW(2)*_I3*(_Cinv(i,j)*_Cinv(k,l) - 0.5*(_Cinv(i,k)*_Cinv(j,l) + _Cinv(i,l)*_Cinv(j,k)) );

    return 4*Cijkl;
  }

} // end namespace GRINS

#endif // GRINS_CARTESIAN_HYPERELASTICITY_H
