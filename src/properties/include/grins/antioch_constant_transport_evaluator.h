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


#ifndef GRINS_ANTIOCH_CONSTANT_TRANSPORT_EVALUATOR_H
#define GRINS_ANTIOCH_CONSTANT_TRANSPORT_EVALUATOR_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// GRINS
#include "grins/antioch_evaluator.h"
#include "grins/antioch_mixture_averaged_transport_mixture.h"

// libMesh
#include "libmesh/libmesh_common.h"

// Antioch
#include "antioch/wilke_evaluator.h"

namespace GRINS
{
  //! Wrapper class for evaluating constant transport properties, including Antioch::ConstantLewisDiffusivity
  /*!
    This class is expected to be constructed *after* threads have been forked and will only
    live during the lifetime of the thread.
    By default, Antioch is working in SI units. Note that this documentation will always
    be built regardless if Antioch is included in the GRINS build or not. Check configure
    output to confirm that Antioch was included in the build.
  */
  template<typename KineticsThermoCurveFit, typename Thermo, typename Conductivity>
  class AntiochConstantTransportEvaluator : public AntiochEvaluator<KineticsThermoCurveFit,Thermo>
  {
  public:

    AntiochConstantTransportEvaluator( const AntiochConstantTransportMixture<KineticsThermoCurveFit,Conductivity>& mixture )
      : AntiochEvaluator<KineticsThermoCurveFit,Thermo>( mixture ),
      _mu( mixture.mu() ),
      _conductivity( mixture.conductivity() ),
      _diffusivity( mixture.diffusivity() )
    {}

    virtual ~AntiochConstantTransportEvaluator() = default;

    libMesh::Real mu( const libMesh::Real /*T*/,
                      const libMesh::Real /*P*/,
                      const std::vector<libMesh::Real>& /*Y*/ )
    { return _mu; }

    libMesh::Real k( const libMesh::Real T,
                     const libMesh::Real /*P*/,
                     const std::vector<libMesh::Real>& Y )
    {
      // Second T is dummy
      const libMesh::Real cp = this->cp( T, T, Y );
      return _conductivity( _mu, cp );
    }

    void mu_and_k_and_D( const libMesh::Real /*T*/,
                         const libMesh::Real rho,
                         const libMesh::Real cp,
                         const std::vector<libMesh::Real>& /*Y*/,
                         libMesh::Real & mu, libMesh::Real & k,
                         std::vector<libMesh::Real> & D )
    {
      mu = _mu;
      k = _conductivity( _mu, cp );
      std::fill( D.begin(), D.end(), _diffusivity.D(rho,cp,k) );
    }

  protected:

    const libMesh::Real _mu;

    const Conductivity & _conductivity;

    const Antioch::ConstantLewisDiffusivity<libMesh::Real> & _diffusivity;

  private:

    AntiochConstantTransportEvaluator();

  };

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_CONSTANT_TRANSPORT_EVALUATOR_H
