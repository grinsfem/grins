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
#include "grins/hookes_law.h"

// GRINS
#include "grins/common.h"
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/tensor_value.h"

namespace GRINS
{
  HookesLaw::HookesLaw(const GetPot& input)
    : StressStrainLaw<HookesLaw>(),
    ParameterUser("HookesLaw"),
    _C(),
    _lambda(0.0),
    _mu(0.0)
  {
    // Warning about this constructor being deprecated
    {
      std::string warning = "WARNING: Use of this constructor is DEPRECATED.\n";
      warning += "         Please update to use constructor with input material name.\n";
      grins_warning(warning);
    }

    this->read_input_options(input);
  }

  HookesLaw::HookesLaw(const GetPot& input, const std::string& material)
    : StressStrainLaw<HookesLaw>(),
    ParameterUser("HookesLaw"),
    _C(),
    _lambda(0.0),
    _mu(0.0)
  {
    MaterialsParsing::duplicate_input_test(input,
                                           "Materials/"+material+"/StressStrainLaw/HookesLaw/lambda",
                                           "Physics/HookesLaw/lambda");
    MaterialsParsing::duplicate_input_test(input,
                                           "Materials/"+material+"/StressStrainLaw/HookesLaw/mu",
                                           "Physics/HookesLaw/mu");
    MaterialsParsing::duplicate_input_test(input,
                                           "Materials/"+material+"/StressStrainLaw/HookesLaw/E",
                                           "Physics/HookesLaw/E");
    MaterialsParsing::duplicate_input_test(input,
                                           "Materials/"+material+"/StressStrainLaw/HookesLaw/nu",
                                           "Physics/HookesLaw/nu");

    // Parse the new version
    if( input.have_variable("Materials/"+material+"/StressStrainLaw/HookesLaw/lambda") &&
        input.have_variable("Materials/"+material+"/StressStrainLaw/HookesLaw/mu") )
      {
        this->set_parameter
          (_lambda, input, "Materials/"+material+"/StressStrainLaw/HookesLaw/lambda", 0.0);
        this->set_parameter
          (_mu, input, "Materials/"+material+"/StressStrainLaw/HookesLaw/mu", 0.0);
      }
    else if( input.have_variable("Materials/"+material+"/StressStrainLaw/HookesLaw/E") &&
             input.have_variable("Materials/"+material+"/StressStrainLaw/HookesLaw/nu") )
      {
        /*! \todo we'll need a special accessor to give ParameterUser access to these */
        libMesh::Real E  = input("Materials/"+material+"/StressStrainLaw/HookesLaw/E", 0.0);
        libMesh::Real nu = input("Materials/"+material+"/StressStrainLaw/HookesLaw/nu", 0.0);
        _lambda = nu*E/( (1+nu)*(1-2*nu) );
        _mu = E/(2*(1+nu));
      }
    // Parse the old version
    else if( input.have_variable("Physics/HookesLaw/lambda") &&
             input.have_variable("Physics/HookesLaw/mu") )
      {
        MaterialsParsing::dep_input_warning( "Physics/HookesLaw/lambda",
                                             "StressStrainLaw/HookesLaw/lambda" );
        MaterialsParsing::dep_input_warning( "Physics/HookesLaw/mu",
                                             "StressStrainLaw/HookesLaw/mu" );

        this->set_parameter
          (_lambda, input, "Physics/HookesLaw/lambda", 0.0);
        this->set_parameter
          (_mu, input, "Physics/HookesLaw/mu", 0.0);
      }
    else if( input.have_variable("Physics/HookesLaw/E") &&
             input.have_variable("Physics/HookesLaw/nu") )
      {
        MaterialsParsing::dep_input_warning( "Physics/HookesLaw/E",
                                             "StressStrainLaw/HookesLaw/E" );
        MaterialsParsing::dep_input_warning( "Physics/HookesLaw/nu",
                                             "StressStrainLaw/HookesLaw/nu" );

        /*! \todo we'll need a special accessor to give ParameterUser access to these */
        libMesh::Real E  = input("Physics/HookesLaw/E", 0.0);
        libMesh::Real nu = input("Physics/HookesLaw/nu", 0.0);
        _lambda = nu*E/( (1+nu)*(1-2*nu) );
        _mu = E/(2*(1+nu));
      }
    else
      {
        libmesh_error_msg("ERROR: Could not find consistent HookesLaw input!");
      }

    // mu should be positive
    if( _mu <= 0.0 )
      libmesh_error_msg("ERROR: Detected non-positive shear modulus!");

    // Technically, lambda can be negative, but it's weird. Let's warn about it.
    if( _lambda <= 0.0 )
      {
        std::string warning = "WARNING: Detected non-positive Lame parameter";
        grins_warning(warning);
      }
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
      this->set_parameter
        (_lambda, input, "Physics/HookesLaw/lambda", 0.0);

    if( input.have_variable("Physics/HookesLaw/mu") )
      this->set_parameter
        (_mu, input, "Physics/HookesLaw/mu", 0.0);

    if( input.have_variable("Physics/HookesLaw/E") &&
        input.have_variable("Physics/HookesLaw/nu") )
      {
        // FIXME - we'll need a special accessor to give parameter
        // access to these
        libMesh::Real E  = input("Physics/HookesLaw/E", 0.0);
        libMesh::Real nu = input("Physics/HookesLaw/nu", 0.0);
        _lambda = nu*E/( (1+nu)*(1-2*nu) );
        _mu = E/(2*(1+nu));
      }
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
  }

  libMesh::Real HookesLaw::compute_33_stress_imp( const libMesh::TensorValue<libMesh::Real>& g_contra,
                                                  const libMesh::TensorValue<libMesh::Real>& g_cov,
                                                  const libMesh::TensorValue<libMesh::Real>& /*G_contra*/,
                                                  const libMesh::TensorValue<libMesh::Real>& G_cov )
  {
    libMesh::Real sigma_33 = 0.0;

    for( unsigned int k = 0; k < 3; k++ )
      {
        for( unsigned int l = 0; l < 3; l++ )
          {
            libMesh::Real strain_kl = 0.5*(G_cov(k,l) - g_cov(k,l));

            libMesh::Real C = _lambda*g_contra(2,2)*g_contra(k,l) +
              _mu*(g_contra(2,k)*g_contra(2,l) + g_contra(2,l)*g_contra(2,k));

            sigma_33 += C*strain_kl;
          }
      }

    return sigma_33;
  }

} // end namespace GRINS
