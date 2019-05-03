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

// This class
#include "grins/cartesian_hyperelasticity.h"

namespace GRINS
{
  template<typename StrainEnergy>
  CartesianHyperlasticity<StrainEnergy>::CartesianHyperlasticity
  ( const libMesh::Tensor & F, const HyperelasticStrainEnergy<StrainEnergy> & W )
    : CartesianStressStrainLaw<CartesianHyperlasticity<StrainEnergy>>(F),
    _W(W),
    _C(this->right_cauchy_green(this->_F)),
    _Cinv(_C.inverse())
  {
    this->compute_invariants(_C,_I1,_I2,_I3);
    _dW(0) = _W.dI1(_I1,_I2,_I3);
    _dW(1) = _W.dI2(_I1,_I2,_I3);
    _dW(2) = _W.dI3(_I1,_I2,_I3);

    // We only need upper triangular part of the second derivative
    _d2W(0,0) = _W.dI12(_I1,_I2,_I3);
    _d2W(1,1) = _W.dI22(_I1,_I2,_I3);
    _d2W(2,2) = _W.dI32(_I1,_I2,_I3);

    _d2W(0,1) = _W.dI1dI2(_I1,_I2,_I3);
    _d2W(0,2) = _W.dI1dI3(_I1,_I2,_I3);
    _d2W(1,2) = _W.dI2dI3(_I1,_I2,_I3);

    _S = this->compute_pk2_stress();
  }

  template<typename StrainEnergy>
  libMesh::Tensor CartesianHyperlasticity<StrainEnergy>::compute_pk2_stress() const
  {
    libMesh::Tensor I(1,0,0,0,1,0,0,0,1);
    return 2*( _dW(0)*I + _dW(1)*(_I1*I - _C) + (_dW(2)*_I3)*_Cinv );
  }

} // end namespace GRINS
