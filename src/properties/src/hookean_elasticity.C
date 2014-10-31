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

// This class
#include "grins/hookean_elasticity.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/tensor_value.h"

namespace GRINS
{
  HookeanElasiticty::HookeanElasiticty(const GetPot& input)
    : ElasticityTensor(),
      _lambda(0.0),
      _mu(0.0)
  {
    this->read_input_options(input);

    // Build initial elasicity tensor using Kronecker delta
    libMesh::TensorValue<libMesh::Real> delta(1.0,0.0,0.0,
                                              0.0,1.0,0.0,
                                              0.0,0.0,1.0);
    this->recompute_elasticity(delta);

    return;
  }

  HookeanElasiticty::~HookeanElasiticty()
  {
    return;
  }

  void HookeanElasiticty::read_input_options(const GetPot& input)
  {
    // We'd better have either Lam\'{e} constants or E and nu
    if( ( !input.have_variable("Physics/HookeanElasticity/lambda") || 
          !input.have_variable("Physics/HookeanElasticity/mu") ) &&
        ( !input.have_variable("Physics/HookeanElasticity/E") || 
          !input.have_variable("Physics/HookeanElasticity/nu") ) )
      {
        std::cerr << "Error: Must specify either Lame constants lambda and mu or" << std::endl
                  << "       Young's modulus and Poisson's ratio." << std::endl;
        libmesh_error();
      }

    if( input.have_variable("Physics/HookeanElasticity/lambda") )
      _lambda = input("Physics/HookeanElasticity/lambda", 0.0 );
    
    if( input.have_variable("Physics/HookeanElasticity/mu") )
      _mu = input("Physics/HookeanElasticity/mu", 0.0 );
        
    if( input.have_variable("Physics/HookeanElasticity/E") && 
        input.have_variable("Physics/HookeanElasticity/nu") )
      {
        libMesh::Real E  = input("Physics/HookeanElasticity/E", 0.0);
        libMesh::Real nu = input("Physics/HookeanElasticity/nu", 0.0);
        _lambda = nu*E/( (1+nu)*(1-2*nu) );
        _mu = E/(2*(1+nu));
      }

    return;
  }
  
  void HookeanElasiticty::recompute_elasticity( libMesh::TensorValue<libMesh::Real>& g )
  {
    for( unsigned int i = 0; i < 3; i++ )
      {
        for( unsigned int j = 0; j < 3; j++ )
          {
            for( unsigned int k = 0; k < 3; k++ )
              {
                for( unsigned int l = 0; l < 3; l++ )
                  {
                    _C[i][j][k][l] = _lambda*g(i,j)*g(k,l) + _mu*(g(i,k)*g(j,l) + g(i,l)*g(j,k));
                  }
              }
          }
      }

    return;
  }

} // end namespace GRINS
