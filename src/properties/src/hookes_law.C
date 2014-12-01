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
#include "grins/hookes_law.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/tensor_value.h"

namespace GRINS
{
  HookesLaw::HookesLaw(const GetPot& input)
    : StressStrainLaw<HookesLaw>(),
      _C(),
      _lambda(0.0),
      _mu(0.0)
  {
    this->read_input_options(input);

    return;
  }

  HookesLaw::~HookesLaw()
  {
    return;
  }

  void HookesLaw::read_input_options(const GetPot& input)
  {
    // We'd better have either Lam\'{e} constants or E and nu
    if( ( !input.have_variable("Physics/HookesLaw/lambda") || 
          !input.have_variable("Physics/HookesLaw/mu") ) &&
        ( !input.have_variable("Physics/HookesLaw/E") || 
          !input.have_variable("Physics/HookesLaw/nu") ) )
      {
        std::cerr << "Error: Must specify either Lame constants lambda and mu or" << std::endl
                  << "       Young's modulus and Poisson's ratio." << std::endl;
        libmesh_error();
      }

    if( input.have_variable("Physics/HookesLaw/lambda") )
      _lambda = input("Physics/HookesLaw/lambda", 0.0 );
    
    if( input.have_variable("Physics/HookesLaw/mu") )
      _mu = input("Physics/HookesLaw/mu", 0.0 );
        
    if( input.have_variable("Physics/HookesLaw/E") && 
        input.have_variable("Physics/HookesLaw/nu") )
      {
        libMesh::Real E  = input("Physics/HookesLaw/E", 0.0);
        libMesh::Real nu = input("Physics/HookesLaw/nu", 0.0);
        _lambda = nu*E/( (1+nu)*(1-2*nu) );
        _mu = E/(2*(1+nu));
      }

    return;
  }
  
  void HookesLaw::compute_stress_imp( unsigned int dim,
                                      const libMesh::TensorValue<libMesh::Real>& g_contra,
                                      const libMesh::TensorValue<libMesh::Real>& g_cov,
                                      const libMesh::TensorValue<libMesh::Real>& /*G_contra*/,
                                      const libMesh::TensorValue<libMesh::Real>& G_cov,
                                      libMesh::TensorValue<libMesh::Real>& stress )
  {
    stress.zero();

    for( unsigned int i = 0; i < dim; i++ )
      {
        for( unsigned int j = 0; j < dim; j++ )
          {
            for( unsigned int k = 0; k < dim; k++ )
              {
                for( unsigned int l = 0; l < dim; l++ )
                  {
                    libMesh::Real strain_kl = 0.5*(G_cov(k,l) - g_cov(k,l));

                    _C(i,j,k,l) = _lambda*g_contra(i,j)*g_contra(k,l) +
                                  _mu*(g_contra(i,k)*g_contra(j,l) + g_contra(i,l)*g_contra(j,k));

                    stress(i,j) += _C(i,j,k,l)*strain_kl;
                  }
              }
          }
      }

    return;
  }

  void HookesLaw::compute_stress_and_elasticity_imp( unsigned int dim,
                                                     const libMesh::TensorValue<libMesh::Real>& g_contra,
                                                     const libMesh::TensorValue<libMesh::Real>& g_cov,
                                                     const libMesh::TensorValue<libMesh::Real>& G_contra,
                                                     const libMesh::TensorValue<libMesh::Real>& G_cov,
                                                     libMesh::TensorValue<libMesh::Real>& stress,
                                                     ElasticityTensor& C)
  {
    this->compute_stress_imp(dim,g_contra,g_cov,G_contra,G_cov,stress);

    C = _C;

    return;
  }

  libMesh::Real HookesLaw::compute_33_stress_imp( const libMesh::TensorValue<libMesh::Real>& g_contra,
                                                  const libMesh::TensorValue<libMesh::Real>& g_cov,
                                                  const libMesh::TensorValue<libMesh::Real>& G_contra,
                                                  const libMesh::TensorValue<libMesh::Real>& G_cov )
  {
    libmesh_not_implemented();
    return 0.0;
  }

} // end namespace GRINS
