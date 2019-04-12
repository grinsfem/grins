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
#include "grins/hookes_law_1d.h"

// GRINS
#include "grins/elasticity_tensor.h"
#include "grins/common.h"
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/tensor_value.h"

namespace GRINS
{
  HookesLaw1D::HookesLaw1D(const GetPot& input)
    : StressStrainLaw<HookesLaw1D>(),
    ParameterUser("HookesLaw1D"),
    _E(0.0),
    _nu(0.0)
  {
    // Warning about this constructor being deprecated
    {
      std::string warning = "WARNING: Use of this constructor is DEPRECATED.\n";
      warning += "         Please update to use constructor with input material name.\n";
      grins_warning(warning);
    }

    this->read_input_options(input);

    return;
  }

  HookesLaw1D::HookesLaw1D(const GetPot& input, const std::string& material)
    : StressStrainLaw<HookesLaw1D>(),
    ParameterUser("HookesLaw1D"),
    _E(0.0),
    _nu(0.0)
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
        // FIXME - we'll need a special accessor to give parameter
        // access to these
        libMesh::Real lambda  = input("Physics/HookesLaw/lambda", 0.0);
        libMesh::Real mu = input("Physics/HookesLaw/mu", 0.0);
        _E  = mu*(3*lambda + 2*mu)/(lambda+mu);
        _nu = lambda/(2*(lambda+mu));
      }
    else if( input.have_variable("Materials/"+material+"/StressStrainLaw/HookesLaw/E") &&
             input.have_variable("Materials/"+material+"/StressStrainLaw/HookesLaw/nu") )
      {
        this->set_parameter
          (_E, input, "Materials/"+material+"/StressStrainLaw/HookesLaw/E", 0.0);
        this->set_parameter
          (_nu, input, "Materials/"+material+"/StressStrainLaw/HookesLaw/nu", 0.0);
      }
    // Parse the old version
    else if( input.have_variable("Physics/HookesLaw/lambda") &&
             input.have_variable("Physics/HookesLaw/mu") )
      {
        MaterialsParsing::dep_input_warning( "Physics/HookesLaw/lambda",
                                             "StressStrainLaw/HookesLaw/lambda" );
        MaterialsParsing::dep_input_warning( "Physics/HookesLaw/mu",
                                             "StressStrainLaw/HookesLaw/mu" );

        /*! \todo we'll need a special accessor to give ParameterUser access to these */
        libMesh::Real lambda  = input("Physics/HookesLaw/lambda", 0.0);
        libMesh::Real mu = input("Physics/HookesLaw/mu", 0.0);
        _E  = mu*(3*lambda + 2*mu)/(lambda+mu);
        _nu = lambda/(2*(lambda+mu));
      }
    else if( input.have_variable("Physics/HookesLaw/E") &&
             input.have_variable("Physics/HookesLaw/nu") )
      {
        MaterialsParsing::dep_input_warning( "Physics/HookesLaw/E",
                                             "StressStrainLaw/HookesLaw/E" );
        MaterialsParsing::dep_input_warning( "Physics/HookesLaw/nu",
                                             "StressStrainLaw/HookesLaw/nu" );

        this->set_parameter(_E, input, "Physics/HookesLaw/E", 0.0);
        this->set_parameter(_nu, input, "Physics/HookesLaw/nu", 0.0);
      }
    else
      {
        libmesh_error_msg("ERROR: Could not find consistent HookesLaw input!");
      }

    // mu should be positive
    if( _E <= 0.0 )
      libmesh_error_msg("ERROR: Detected non-positive Young's modulus!");

    // Technically, Poisson's ratio can be negative, but it's weird. Let's warn about it.
    if( _nu < 0.0 )
      {
        std::string warning = "WARNING: Detected non-positive Poisson's ratio";
        grins_warning(warning);
      }
  }

  HookesLaw1D::~HookesLaw1D()
  {
    return;
  }

  void HookesLaw1D::read_input_options(const GetPot& input)
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

    if( input.have_variable("Physics/HookesLaw/lambda") &&
        input.have_variable("Physics/HookesLaw/mu") )
      {
        // FIXME - we'll need a special accessor to give parameter
        // access to these
        libMesh::Real lambda  = input("Physics/HookesLaw/lambda", 0.0);
        libMesh::Real mu = input("Physics/HookesLaw/mu", 0.0);
        _E  = mu*(3*lambda + 2*mu)/(lambda+mu);
        _nu = lambda/(2*(lambda+mu));
      }
    else
      {
        if( input.have_variable("Physics/HookesLaw/E") )
          this->set_parameter
            (_E, input, "Physics/HookesLaw/E", 0.0);

        if( input.have_variable("Physics/HookesLaw/nu") )
          this->set_parameter
            (_nu, input, "Physics/HookesLaw/nu", 0.0);
      }
  }

  void HookesLaw1D::compute_stress_imp( unsigned int /*dim*/,
                                        const libMesh::TensorValue<libMesh::Real>& g_contra,
                                        const libMesh::TensorValue<libMesh::Real>& g_cov,
                                        const libMesh::TensorValue<libMesh::Real>& /*G_contra*/,
                                        const libMesh::TensorValue<libMesh::Real>& G_cov,
                                        libMesh::TensorValue<libMesh::Real>& stress )
  {
    stress.zero();

    libMesh::Real strain = 0.5*(G_cov(0,0) - g_cov(0,0));

    stress(0,0) = (this->_E)*g_contra(0,0)*g_contra(0,0)*strain;

    return;
  }

  void HookesLaw1D::compute_stress_and_elasticity_imp( unsigned int dim,
                                                       const libMesh::TensorValue<libMesh::Real>& g_contra,
                                                       const libMesh::TensorValue<libMesh::Real>& g_cov,
                                                       const libMesh::TensorValue<libMesh::Real>& G_contra,
                                                       const libMesh::TensorValue<libMesh::Real>& G_cov,
                                                       libMesh::TensorValue<libMesh::Real>& stress,
                                                       ElasticityTensor& C)
  {
    this->compute_stress_imp(dim,g_contra,g_cov,G_contra,G_cov,stress);

    C(0,0,0,0) = this->_E*g_contra(0,0)*g_contra(0,0);

    return;
  }

  libMesh::Real HookesLaw1D::compute_33_stress_imp( const libMesh::TensorValue<libMesh::Real>& /*g_contra*/,
                                                    const libMesh::TensorValue<libMesh::Real>& /*g_cov*/,
                                                    const libMesh::TensorValue<libMesh::Real>& /*G_contra*/,
                                                    const libMesh::TensorValue<libMesh::Real>& /*G_cov*/ )
  {
    libmesh_not_implemented();

    return 0.0;
  }

} // end namespace GRINS
