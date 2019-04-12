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


#ifndef GRINS_ANTIOCH_WILKE_TRANSPORT_EVALUATOR_H
#define GRINS_ANTIOCH_WILKE_TRANSPORT_EVALUATOR_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// GRINS
#include "grins/antioch_evaluator.h"
#include "grins/antioch_mixture_averaged_transport_mixture.h"

// libMesh
#include "libmesh/libmesh_common.h"

// Antioch
#include "antioch/mixture_averaged_transport_evaluator.h"

namespace GRINS
{
  //! Wrapper class for evaluating Wilke transport properties using Antioch
  /*!
    This class is expected to be constructed *after* threads have been forked and will only
    live during the lifetime of the thread.
    By default, Antioch is working in SI units. Note that this documentation will always
    be built regardless if Antioch is included in the GRINS build or not. Check configure
    output to confirm that Antioch was included in the build.
  */
  template<typename KineticsThermoCurveFit, typename Thermo, typename Viscosity, typename Conductivity, typename Diffusivity>
  class AntiochMixtureAveragedTransportEvaluator : public AntiochEvaluator<KineticsThermoCurveFit,Thermo>
  {
  public:

    AntiochMixtureAveragedTransportEvaluator( const AntiochMixtureAveragedTransportMixture<KineticsThermoCurveFit,Thermo,Viscosity,Conductivity,Diffusivity>& mixture );

    virtual ~AntiochMixtureAveragedTransportEvaluator(){};

    libMesh::Real mu( const libMesh::Real T,
                      const libMesh::Real P,
                      const std::vector<libMesh::Real>& Y );

    libMesh::Real k( const libMesh::Real T,
                     const libMesh::Real P,
                     const std::vector<libMesh::Real>& Y );

    void mu_and_k_and_D( const libMesh::Real T,
                         const libMesh::Real rho,
                         const libMesh::Real cp,
                         const std::vector<libMesh::Real>& Y,
                         libMesh::Real& mu, libMesh::Real& k,
                         std::vector<libMesh::Real>& D );


  protected:

    std::unique_ptr<Antioch::MixtureAveragedTransportEvaluator<Diffusivity,Viscosity,Conductivity,libMesh::Real> > _wilke_evaluator;

    const Antioch::MixtureDiffusion<Diffusivity,libMesh::Real>& _diffusivity;

  private:

    AntiochMixtureAveragedTransportEvaluator();

  };

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_WILKE_TRANSPORT_EVALUATOR_H
