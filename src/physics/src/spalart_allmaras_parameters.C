//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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
#include "grins/spalart_allmaras_parameters.h"

namespace GRINS
{
  SpalartAllmarasParameters::SpalartAllmarasParameters(const GetPot& input )
    : _kappa(0.41),
      _cv1(7.1),
      _cv2(0.7),
      _cv3(0.9),
      _cb1(0.1355),
      _sigma(2./3.),
      _cb2(0.622),
      _cw1( _cb1/(_kappa*_kappa) + (1.0+_cb2)/_sigma ),
      _r_lin(10.0),
      _c_w2(0.3),
      _c_w3(2.0),
      _c_t3(1.2),
      _c_t4(0.5),
      _c_n1(16.0)
  {}


  libMesh::Real SpalartAllmarasParameters::source_fn(libMesh::Number nu, libMesh::Real mu,
                                                     libMesh::Real wall_distance, libMesh::Real vorticity_value) const
  {
    // Step 1
    libMesh::Real chi = nu/mu;

    // Step 2
    libMesh::Real chi3 = chi*chi*chi;
    libMesh::Real fv1 = chi3/(chi3 + _cv1*_cv1*_cv1);

    // Step 3
    libMesh::Real fv2 = 1 - (chi/(1 + chi*fv1));

    // Step 4
    libMesh::Real S_bar = nu/(pow(_kappa, 2.0) * pow(wall_distance, 2.0))*(fv2) ;

    // Step 5, the absolute value of the vorticity
    libMesh::Real S = vorticity_value;

    // Step 6
    libMesh::Real S_tilde = 0.0;
    if(S_bar >= -this->_cv2*S)
      {
        S_tilde = S + S_bar;
      }
    else
      {
        S_tilde = S + (S*(pow(this->_cv2,2.0)*S + this->_cv3*S_bar))/((this->_cv3 - (2*this->_cv2))*S - S_bar);
      }

    return S_tilde;
  }

  libMesh::Real SpalartAllmarasParameters::destruction_fn(libMesh::Number nu, libMesh::Real wall_distance,
                                                          libMesh::Real S_tilde) const
  {
    // Step 1
    libMesh::Real r = 0.0;

    r = std::min(nu/(S_tilde*pow(this->_kappa,2.0)*pow(wall_distance,2.0)), this->_r_lin);

    // Step 2
    libMesh::Real g = r + this->_c_w2*(pow(r,6.0) - r);

    // Step 3
    libMesh::Real fw = g*pow((1 + pow(this->_c_w3,6.0))/(pow(g,6.0) + pow(this->_c_w3,6.0)), 1.0/6.0);

    return fw;
  }

} // end namespace GRINS
