//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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

// GRINS
#include "grins/physics_naming.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  SpalartAllmarasParameters::SpalartAllmarasParameters(const GetPot& input )
    : ParameterUser("SpalartAllmarasParameters"),
      _kappa(input("Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/kappa",0.41)),
      _cv1(input("Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/cv1",7.1)),
      _cv2(input("Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/cv2",0.7)),
      _cv3(input("Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/cv3",0.9)),
      _cb1(input("Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/cb1",0.1355)),
      _sigma(input("Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/sigma",2./3.)),
      _cb2(input("Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/cb2",0.622)),
      //_cw1( _cb1/(_kappa*_kappa) + (1.0+_cb2)/_sigma ),
      _r_lin(input("Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/r_lin",10.0)),
      _c_w2(input("Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/c_w2",0.3)),
      _c_w3(input("Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/c_w3",2.0)),
      _c_t3(input("Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/c_t3",1.2)),
      _c_t4(input("Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/c_t4",0.5)),
      _c_n1(input("Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/c_n1",16.0))
  {
    this->set_parameter(this->_kappa ,input, "Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/kappa" , this->_kappa );
    this->set_parameter(this->_cv1 ,input, "Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/cv1", this->_cv1 );
    this->set_parameter(this->_cv2 ,input, "Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/cv2", this->_cv2 );
    this->set_parameter(this->_cv3 ,input, "Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/cv3", this->_cv3 );
    this->set_parameter(this->_cb1 ,input, "Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/cb1", this->_cb1 );
    this->set_parameter(this->_sigma ,input, "Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/sigma", this->_sigma );
    this->set_parameter(this->_cb2 ,input, "Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/cb2" , this->_cb2 );
    this->set_parameter(this->_r_lin ,input, "Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/r_lin", this->_r_lin );
    this->set_parameter(this->_c_w2 ,input, "Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/c_w2", this->_c_w2 );
    this->set_parameter(this->_c_w3 ,input, "Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/c_w3", this->_c_w3 );
    this->set_parameter(this->_c_t3 ,input, "Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/c_t3", this->_c_t3 );
    this->set_parameter(this->_c_t4 ,input, "Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/c_t4", this->_c_t4 );
    this->set_parameter(this->_c_n1 ,input, "Physics/"+PhysicsNaming::spalart_allmaras()+"/Parameters/c_n1", this->_c_n1 );
  }

  libMesh::Real SpalartAllmarasParameters::source_fn(libMesh::Number nu, libMesh::Real mu,
                                                     libMesh::Real wall_distance, libMesh::Real vorticity_value, bool infinite_distance) const
  {
    // Step 1
    libMesh::Real chi = nu/mu;

    // Step 2
    libMesh::Real fv1 = this->fv1(chi);

    // Step 3
    libMesh::Real fv2 = 1 - (chi/(1 + chi*fv1));

    libMesh::Real S_bar = 0.0;
    // Step 4
    if(infinite_distance)
      {
        S_bar = 0.0;
      }
    else
      {
        S_bar = nu/(pow(_kappa, 2.0) * pow(wall_distance, 2.0))*(fv2) ;
      }

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
                                                          libMesh::Real S_tilde, bool infinite_distance) const
  {
    // Step 1
    libMesh::Real r = 0.0;
    if(infinite_distance)
      {
        r = 0.0;
      }
    else
      {
        r = std::min(nu/(S_tilde*pow(this->_kappa,2.0)*pow(wall_distance,2.0)), this->_r_lin);
      }

    // Step 2
    libMesh::Real g = r + this->_c_w2*(pow(r,6.0) - r);

    // Step 3
    libMesh::Real fw = g*pow(((1 + pow(this->_c_w3,6.0))/(pow(g,6.0) + pow(this->_c_w3,6.0))), 1.0/6.0);

    return fw;
  }

} // end namespace GRINS
